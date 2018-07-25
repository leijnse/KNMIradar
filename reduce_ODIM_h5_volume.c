/*This program reduces the radar volume files by: */
/* 1) conversion from 16 to 8-bit values          */
/* 2) reducing the number of moments              */
/* 3) reducing the number of scans                */
/* 4) reducing the range resolution               */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <zlib.h>

herr_t copy_attributes(hid_t loc_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data);
int resample_range(unsigned short *data16, long nbins, long nrays, double nodata, double undetect, double rscale, double res, long nbins_new, unsigned short *data16_new);
int convert_16_to_8bit(unsigned short *data16, long nrays, long nbins, double *gain, double *offset, double *nodata, double *undetect, unsigned char *data8);
hid_t H5LTmake_dataset_zip(hid_t loc_id, char *dset_name, int rank, hsize_t *dims, hid_t tid, void *data);

#define LSTR 256
#define NMOMENTX 256

int main(int argc, char **argv)
{
	int i, j, i_data, i_data_new, Nscan, scan, scan_new, Nrmmom, remove_moment, *keep_scan, *to_8bit;
	long nbins, nrays, nbins_new;
	double gain, offset, nodata, undetect, tmp_res, rscale, *resample_res;
	char object[LSTR], object_new[LSTR], string[LSTR], rmmom[NMOMENTX][LSTR];
	unsigned short *data16, *data16_new;
	unsigned char *data8;
	hid_t h5_in, h5_out, gid, did;
	hsize_t dims[2];
	H5T_class_t class_id;
	size_t type_size;
	
	/*Check if number of input arguments is sufficient.*/
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <infile.h5> <outfile.h5> [-rms<n>/-rmd<S>/-resample<n,x>/-resample<x>/-8bit<n>/-8bit]\n", argv[0]);
		return 1;
	}
	
	/*Open HDF5 input file.*/
	if ((h5_in = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
	{
		fprintf(stderr, "Error opening file %s\n", argv[1]);
		return 2;
	}
	
	/*Count number of dataset groups.*/
	sprintf(object, "/dataset1");
	Nscan = 0;
	while (H5Lexists(h5_in, object, H5P_DEFAULT))
	{
		Nscan += 1;
		sprintf(object, "/dataset%d", Nscan + 1);
	}
	if (Nscan == 0)
	{
		fprintf(stderr, "Error: file does not contain any dataset groups!\n");
		H5Fclose(h5_in);
		return 3;
	}
	
	/*Open HDF5 input file.*/
	if ((h5_out = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	{
		fprintf(stderr, "Error creating file %s\n", argv[2]);
		H5Fclose(h5_in);
		return 4;
	}
	
	/*Copy all root attributes and groups other than the scan groups.*/
	if (H5Lexists(h5_in, "/what", H5P_DEFAULT)) H5Ocopy(h5_in, "/what", h5_out, "/what", H5P_DEFAULT, H5P_DEFAULT);
	if (H5Lexists(h5_in, "/where", H5P_DEFAULT)) H5Ocopy(h5_in, "/where", h5_out, "/where", H5P_DEFAULT, H5P_DEFAULT);
	if (H5Lexists(h5_in, "/how", H5P_DEFAULT)) H5Ocopy(h5_in, "/how", h5_out, "/how", H5P_DEFAULT, H5P_DEFAULT);
	H5Aiterate_by_name(h5_in, ".", H5_INDEX_NAME, H5_ITER_NATIVE, NULL, copy_attributes, &h5_out, H5P_DEFAULT);
	
	/*Check if scans need to be removed.*/
	keep_scan = (int *) malloc(Nscan * sizeof(int));
	resample_res = (double *) malloc(Nscan * sizeof(double));
	to_8bit = (int *) malloc(Nscan * sizeof(int));
	for (i = 0; i < Nscan; i++)
	{
		keep_scan[i] = 1;
		resample_res[i] = -1.0;
		to_8bit[i] = 0;
	}
	Nrmmom = 0;
	for (i = 3; i < argc; i++)
	{
		if (sscanf(argv[i], "-rms%d", &scan) > 0) keep_scan[scan - 1] = 0;
		if (sscanf(argv[i], "-rmd%s", rmmom[Nrmmom]) > 0) Nrmmom += 1;
		if (strncmp(argv[i], "-resample", 9) == 0)
		{
			if (sscanf(argv[i], "-resample%d,%lf", &scan, &tmp_res) == 2) resample_res[scan - 1] = tmp_res;
			else if (sscanf(argv[i], "-resample%lf", &tmp_res) == 1) for (j = 0; j < Nscan; j++) if (resample_res[j] < 0.0) resample_res[j] = tmp_res;
		}
		if (strncmp(argv[i], "-8bit", 5) == 0)
		{
			if (sscanf(argv[i], "-8bit%d", &scan) == 1) to_8bit[scan - 1] = 1;
			else for (j = 0; j < Nscan; j++) to_8bit[j] = 1;
		}
	}
	
	/*Loop over all scans.*/
	scan_new = 0;
	for (scan = 1; scan <= Nscan; scan++)
	{
		/*Check if scan should be included.*/
		if (keep_scan[scan - 1] == 1)
		{
			/*Set scan objects, and create apropriate group in output file.*/
			scan_new += 1;
			sprintf(object, "/dataset%d", scan);
			sprintf(object_new, "/dataset%d", scan_new);
			gid = H5Gcreate(h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			H5Gclose(gid);
			
			/*Copy dataset attributes to the new group.*/
			sprintf(object, "/dataset%d/what", scan);
			sprintf(object_new, "/dataset%d/what", scan_new);
			if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
			sprintf(object, "/dataset%d/where", scan);
			sprintf(object_new, "/dataset%d/where", scan_new);
			if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
			sprintf(object, "/dataset%d/how", scan);
			sprintf(object_new, "/dataset%d/how", scan_new);
			if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
			
			/*Check if reduction is necessary, and if so, read appropriate metadata and allocate memory.*/
			if ((resample_res[scan - 1] > 0.0) || (to_8bit[scan - 1] == 1))
			{
				sprintf(object, "/dataset%d/where", scan);
				H5LTget_attribute_long(h5_in, object, "nbins", &nbins);
				H5LTget_attribute_long(h5_in, object, "nrays", &nrays);
				H5LTget_attribute_double(h5_in, object, "rscale", &rscale);
				data16 = (unsigned short *) malloc(nbins * nrays * sizeof(unsigned short));
				
				/*Compute new number of range bins if range resampling is desired, and allocate memory for resampled dataset.*/
				if (resample_res[scan - 1] > 0.0)
				{
					nbins_new = (long) ((nbins * rscale) / resample_res[scan - 1] + 0.5);
					data16_new = (unsigned short *) malloc(nbins_new * nrays * sizeof(unsigned short));
					
					/*Adjust nbins and rscale metadata in output file.*/
					sprintf(object_new, "/dataset%d/where", scan_new);
					H5LTset_attribute_long(h5_out, object_new, "nbins", &nbins_new, 1);
					H5LTset_attribute_double(h5_out, object_new, "rscale", &(resample_res[scan - 1]), 1);
				}
				else nbins_new = nbins;
				
				/*Allocate 8-bit dataset.*/
				data8 = (unsigned char *) malloc(nbins_new * nrays * sizeof(unsigned char));
			}
			
			/*Loop over all data groups.*/
			i_data = 1;
			i_data_new = 0;
			sprintf(object, "/dataset%d/data%d", scan, i_data);
			while (H5Lexists(h5_in, object, H5P_DEFAULT))
			{
				/*Check if this moment needs to be removed.*/
				remove_moment = 0;
				if (Nrmmom > 0)
				{
					sprintf(object, "/dataset%d/data%d/what", scan, i_data);
					H5LTget_attribute_string(h5_in, object, "quantity", string);
					for (i = 0; i < Nrmmom; i++) if (strcmp(string, rmmom[i]) == 0) break;
					if (i < Nrmmom) remove_moment = 1;
				}
				
				/*Process moment if desired.*/
				if (remove_moment == 0)
				{
					/*Create new data group.*/
					i_data_new += 1;
					sprintf(object_new, "/dataset%d/data%d", scan_new, i_data_new);
					gid = H5Gcreate(h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					H5Gclose(gid);
					
					/*Copy metadata.*/
					sprintf(object, "/dataset%d/data%d/what", scan, i_data);
					sprintf(object_new, "/dataset%d/data%d/what", scan_new, i_data_new);
					if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
					sprintf(object, "/dataset%d/data%d/where", scan, i_data);
					sprintf(object_new, "/dataset%d/data%d/where", scan_new, i_data_new);
					if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
					sprintf(object, "/dataset%d/data%d/how", scan, i_data);
					sprintf(object_new, "/dataset%d/data%d/how", scan_new, i_data_new);
					if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
					
					/*Read number of bytes in data.*/
					sprintf(object, "/dataset%d/data%d/data", scan, i_data);
					H5LTget_dataset_info(h5_in, object, dims, &class_id, &type_size);
					
					/*Check if reduction is necessary. Otherwise copy.*/
					if ((resample_res[scan - 1] < 0.0) && ((to_8bit[scan - 1] == 0) || (type_size < 2)))
					{
						/*Copy data.*/
						sprintf(object, "/dataset%d/data%d/data", scan, i_data);
						sprintf(object_new, "/dataset%d/data%d/data", scan_new, i_data_new);
						if (H5Lexists(h5_in, object, H5P_DEFAULT)) H5Ocopy(h5_in, object, h5_out, object_new, H5P_DEFAULT, H5P_DEFAULT);
					}
					else
					{
						/*Read data and metadata.*/
						sprintf(object, "/dataset%d/data%d/data", scan, i_data);
						H5LTread_dataset(h5_in, object, H5T_NATIVE_USHORT, data16);
						
						/*Read gain, offset, nodata, and undetect.*/
						sprintf(object, "/dataset%d/data%d/what", scan, i_data);
						H5LTget_attribute_double(h5_in, object, "gain", &gain);
						H5LTget_attribute_double(h5_in, object, "offset", &offset);
						H5LTget_attribute_double(h5_in, object, "nodata", &nodata);
						H5LTget_attribute_double(h5_in, object, "undetect", &undetect);
						
						/*Resample in range if desired.*/
						if (resample_res[scan - 1] > 0.0) resample_range(data16, nbins, nrays, nodata, undetect, rscale, resample_res[scan - 1], nbins_new, data16_new);
						
						/*Convert to 8-bit values if desired.*/
						if ((to_8bit[scan - 1] == 1) && (type_size > 1))
						{
							if (resample_res[scan - 1] > 0.0) convert_16_to_8bit(data16_new, nbins_new, nrays, &gain, &offset, &nodata, &undetect, data8);
							else convert_16_to_8bit(data16, nbins, nrays, &gain, &offset, &nodata, &undetect, data8);
							
							/*Adjust gain, offset, nodata, and undetect metadata in output file.*/
							sprintf(object_new, "/dataset%d/data%d/what", scan_new, i_data_new);
							H5LTset_attribute_double(h5_out, object_new, "gain", &gain, 1);
							H5LTset_attribute_double(h5_out, object_new, "offset", &offset, 1);
							H5LTset_attribute_double(h5_out, object_new, "nodata", &nodata, 1);
							H5LTset_attribute_double(h5_out, object_new, "undetect", &undetect, 1);
						}
						
						/*Copy to 8-bit if original was 8-bit.*/
						if (type_size < 2) for (i = 0; i < (nbins_new * nrays); i++) data8[i] = (unsigned char) data16_new[i];
						
						/*Write new data to file.*/
						dims[0] = nrays;
						dims[1] = nbins_new;
						sprintf(object, "/dataset%d/data%d/data", scan, i_data);
						sprintf(object_new, "/dataset%d/data%d/data", scan_new, i_data_new);
						if ((to_8bit[scan - 1] == 1) || (type_size < 2)) did = H5LTmake_dataset_zip(h5_out, object_new, 2, dims, H5T_NATIVE_UCHAR, data8);
						else did = H5LTmake_dataset_zip(h5_out, object_new, 2, dims, H5T_NATIVE_USHORT, data16_new);
						
						/*Copy data metadata.*/
						H5Aiterate_by_name(h5_in, object, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, copy_attributes, &did, H5P_DEFAULT);
						
						/*Close dataset id.*/
						H5Dclose(did);
					}
				}
				
				/*Next variable.*/
				i_data += 1;
				sprintf(object, "/dataset%d/data%d", scan, i_data);
			}
			
			if ((resample_res[scan - 1] > 0.0) || (to_8bit[scan - 1] == 1))
			{
				free(data8);
				free(data16);
				if (resample_res[scan - 1] > 0.0) free(data16_new);
			}
		}
	}
	
	/*Close files.*/
	H5Fclose(h5_in);
	H5Fclose(h5_out);
	
	/*End of program.*/
	return 0;
}



int resample_range(unsigned short *data16, long nbins, long nrays, double nodata, double undetect, double rscale, double res, long nbins_new, unsigned short *data16_new)
{
	int i, j, j0, j1, k, nd, ud;
	double r0, r1, w, w_x, w_nd, w_ud, x;
	
	/*Create integer nodata and undetect.*/
	nd = (int) nodata;
	ud = (int) undetect;
	
	/*Loop over all azimuths.*/
	for (i = 0; i < nrays; i++)
	{
		r1 = 0.0;
		for (j = 0; j < nbins_new; j++)
		{
			r0 = r1;
			r1 = r0 + res;
			j0 = (int) (r0 / rscale);
			j1 = (int) (r1 / rscale);
			w_nd = 0.0;
			w_ud = 0.0;
			w_x = 0.0;
			x = 0.0;
			for (k = j0; k <= j1; k++)
			{
				w = 1.0;
				if (k == j0) w -= (r0 - k * rscale) / rscale;
				if (k == j1) w -= ((k + 1) * rscale - r1) / rscale;
				
				if (data16[i * nbins + k] == nd) w_nd += w;
				else if (data16[i * nbins + k] == ud) w_ud += w;
				else
				{
					w_x += w;
					x += w * data16[i * nbins + k];
				}
			}
			
			/*Compute average value and map to unsigned short.*/
			if ((w_ud > w_x) && (w_ud >= w_nd)) data16_new[i * nbins_new + j] = (unsigned short) undetect;
			else if ((w_nd > w_x) && (w_nd > w_ud)) data16_new[i * nbins_new + j] = (unsigned short) nodata;
			else
			{
				x /= w_x;
				if (x <= 1.0) data16_new[i * nbins_new + j] = 1;
				else if (x >= 65534.0) data16_new[i * nbins_new + j] = 65534;
				else data16_new[i * nbins_new + j] = (unsigned short) (x + 0.5);
			}
		}
	}
	
	/*End of function.*/
	return 1;
}



int convert_16_to_8bit(unsigned short *data16, long nrays, long nbins, double *gain, double *offset, double *nodata, double *undetect, unsigned char *data8)
{
	double factor, gain_new, offset_new;
	int i, nd, ndn, ud, udn;
	
	/*Compute new values for gain, and offset.*/
	factor = 254.0  / 65534.0;
	gain_new = (*gain) / factor;
	offset_new = (*offset) + (*gain) - gain_new;
	
	/*Set new nodata and undetect values, and create integer nodata and undetect.*/
	ndn = (int) (factor * (*nodata - 1.0) + 1.0);
	udn = (int) (factor * (*undetect - 1.0) + 1.0);
	nd = (int) *nodata;
	ud = (int) *undetect;
	
	/*Loop over all ranges and azimuths.*/
	for (i = 0; i < (nrays * nbins); i++)
	{
		/*Check for nodata and undetect, and set values of the 8-bit data.*/
		if (data16[i] == nd) data8[i] = ndn;
		else if (data16[i] == ud) data8[i] = udn;
		else data8[i] = (unsigned char) (factor * (data16[i] - 1) + 1.5);
	}
	
	/*Set new offset, gain, nodata, and undetect values.*/
	(*offset) = offset_new;
	(*gain) = gain_new;
	(*nodata) = (double) ndn;
	(*undetect) = (double) udn;
	
	/*Exit.*/
	return 0;
}



herr_t copy_attributes(hid_t loc_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data)
{
	hid_t *h5_out;
	hid_t aid, atype, sid;
	int rank, size, i;
	void *att_data;
	hsize_t *dims, *maxdims;
	
	h5_out = op_data;
	
	/*Read attribute and its characteristics.*/
	aid = H5Aopen_by_name(loc_id, ".", attr_name, H5P_DEFAULT, H5P_DEFAULT);
	atype = H5Aget_type(aid);
	sid = H5Aget_space(aid);
	rank = H5Sget_simple_extent_ndims(sid);
	dims = (hsize_t *) malloc(rank * sizeof(hsize_t));
	maxdims = (hsize_t *) malloc(rank * sizeof(hsize_t));
	H5Sget_simple_extent_dims(sid, dims, maxdims);
	size = (int) H5Tget_size(atype);
	for (i = 0; i < rank; i++) size *= dims[i];
	att_data = malloc(size);
	H5Aread(aid, atype, att_data);
	H5Aclose(aid);
	H5Sclose(sid);
	
	/*Write attribute.*/
	sid = H5Screate(H5S_SIMPLE);
	H5Sset_extent_simple(sid, rank, dims, maxdims);
	aid = H5Acreate(*h5_out, attr_name, atype, sid, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(aid, atype, att_data);
	H5Aclose(aid);
	H5Sclose(sid);
	
	/*Free memory.*/
	free(dims);
	free(maxdims);
	free(att_data);
	
	return 0;
}



hid_t H5LTmake_dataset_zip(hid_t loc_id, char *dset_name, int rank, hsize_t *dims, hid_t tid, void *data) 
{
	#define H5ZIPLEVEL (6)
	hid_t sid, pid, did;
	
	/*Create the data space for the dataset.*/
	if ((sid = H5Screate_simple(rank, dims, NULL)) < 0) return -1;
	
	/*Create the property list for zipped datasets.*/
	pid = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(pid, rank, dims);
	H5Pset_deflate(pid, H5ZIPLEVEL);
	
	/*Create the dataset.*/
	if ((did = H5Dcreate(loc_id, dset_name, tid, sid, H5P_DEFAULT, pid, H5P_DEFAULT)) < 0) goto out;
	
	/*Write the dataset only if there is data to write*/
	if (data)
	{
		if (H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) goto out;
	}
	
	/*Terminate access to the data space. */
	if (H5Sclose(sid) < 0) return -1;
	
	/*End access to the property list.*/
	if (H5Pclose(pid) < 0) return -1;
	return did;
	
	out:
		H5Dclose(did);
		H5Sclose(sid);
		H5Sclose(pid);
		return -1;
}
