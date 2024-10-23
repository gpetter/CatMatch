from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.table import Table, vstack, hstack, join
import logging


def get_coordinate_keys(catalog):
	column_names = catalog.colnames

	ra_keys = [s for s in column_names if "ra" in s.lower()]
	dec_keys = [s for s in column_names if "de" in s.lower()]
	if min(len(ra_keys), len(dec_keys)) == 0:
		raise ValueError("No *RA*/*ra* or *DE*/*de* pairs found.")
	ra_key, dec_key = ra_keys[0], dec_keys[0]

	logging.warning(f"'{ra_key}' will be used as RA coordinate")
	logging.warning(f"'{dec_key}' will be used as DEC coordinate")
	return ra_key, dec_key

def match_coords(coords1, coords2, max_sep=1.0, match_type="best"):
	# add an astropy unit of arcsec to matching radius
	max_sep *= u.arcsec
	lons1, lats1 = np.array(coords1[0]), np.array(coords1[1])
	lons2, lats2 = np.array(coords2[0]), np.array(coords2[1])
	# keep track of original indices
	original_idxs1 = np.arange(len(lons1))
	original_idxs2 = np.arange(len(lons2))

	skycoord1 = SkyCoord(np.array(lons1), np.array(lats1), unit='deg', frame='icrs')
	skycoord2 = SkyCoord(np.array(lons2), np.array(lats2), unit='deg', frame='icrs')

	if match_type == "all":
		idx1, idx2, d2d, d3d = skycoord2.search_around_sky(skycoord1, seplimit=max_sep)

		separations = d2d.to("arcsec").value

	else:
		idx2, d2d, d3d = skycoord1.match_to_catalog_sky(skycoord2)
		sep_constraint = np.where(d2d < max_sep)[0]
		idx1 = original_idxs1[sep_constraint]
		idx2 = original_idxs2[idx2[sep_constraint]]

		separations = d2d[sep_constraint].to("arcsec").value

		if match_type == "best_symmetric":
			newidx1, d2d_2, d3d_2 = skycoord2.match_to_catalog_sky(skycoord1)
			sep_constraint = np.where(d2d_2 < max_sep)[0]
			newidx2 = original_idxs2[sep_constraint]
			newidx1 = original_idxs1[newidx1[sep_constraint]]
			new_pairs = np.transpose([newidx1, newidx2])
			old_pairs = np.transpose([idx1, idx2])
			common_pairs = old_pairs[(old_pairs[:, None] == new_pairs).all(-1).any(1)].transpose()
			idx1, idx2 = common_pairs[0], common_pairs[1]

			separations = d2d[idx1].to("arcsec").value
	return idx1, idx2, separations


def match_catalogs(
	catalog_or_coords1,
	catalog_or_coords2,
	max_sep=1.0,
	match_type="best", # ["best", "best_symmetric", "all"]
	join_type="1+2",
	catalog1_coord_names=None,
	catalog2_coord_names=None,
):
	"""
	find matches between two catalogs of coordinates
	Parameters
	----------
	catalog_or_coords1: tuple of (lon, lat) or an astropy Table containing columns 'RA', 'DEC'
	catalog_or_coords2: same as 1
	max_sep: float in units of arcsec
	match_type: str, one of ['best', 'best_symmetric', 'all']

	Returns
	-------

	"""

	if catalog1_coord_names is None:
		catalog1_coord_names = get_coordinate_keys(catalog_or_coords1)
	if catalog2_coord_names is None:
		catalog2_coord_names = get_coordinate_keys(catalog_or_coords2)
	lons1, lats1 = (catalog_or_coords1[catalog1_coord_names[0]],
					catalog_or_coords1[catalog1_coord_names[1]])
	lons2, lats2 = (catalog_or_coords2[catalog2_coord_names[0]],
					catalog_or_coords2[catalog2_coord_names[1]])

	idx1, idx2, separations = match_coords(
		coords1=(lons1, lats1),
		coords2=(lons2, lats2),
		max_sep=max_sep,
		match_type=match_type,
	)


	output_table = hstack((catalog_or_coords1[idx1], catalog_or_coords2[idx2]))
	if match_type == "best":

		unique, counts = np.unique(idx2, return_counts=True)
		counts_per_entry = counts[np.searchsorted(unique, idx2)]
		duplicates = np.where(counts_per_entry > 1)[0]

		duplicate_counts = np.full(len(idx2), np.nan)
		duplicate_counts[duplicates] = counts_per_entry[duplicates]
		output_table["GroupSize"] = duplicate_counts
	if match_type == "all":
		unique1, counts1 = np.unique(idx1, return_counts=True)
		counts_per_entry1 = counts1[np.searchsorted(unique1, idx1)]
		duplicates1 = np.where(counts_per_entry1 > 1)[0]

		unique2, counts2 = np.unique(idx2, return_counts=True)
		counts_per_entry2 = counts2[np.searchsorted(unique2, idx2)]
		duplicates2 = np.where(counts_per_entry2 > 1)[0]

		duplicate_counts = np.full(len(idx2), np.nan)
		duplicate_counts[duplicates1] = counts_per_entry1[duplicates1]
		duplicate_counts[duplicates2] = counts_per_entry2[duplicates2]
		output_table["GroupSize"] = duplicate_counts

	output_table["Separation"] = separations


	if join_type == "1+2":

		return output_table
	elif join_type == "all_from_1":
		nomatch_1 = catalog_or_coords1[np.setdiff1d(np.arange(len(catalog_or_coords1)), idx1)]
		return vstack((output_table, nomatch_1))
	elif join_type == "1_not_2":
		nomatch_1 = catalog_or_coords1[np.setdiff1d(np.arange(len(catalog_or_coords1)), idx1)]
		return nomatch_1[catalog_or_coords1.colnames]



