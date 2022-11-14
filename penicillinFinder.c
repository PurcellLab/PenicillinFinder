#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXSTR 50000
#define MAXIONS 2000
#define MAXPROT 100000
#define FRAGTOL 0.05
#define PRECTOL 0.05

typedef struct {
	int cnt;
	int *prot_ids;
} peptide_evidence;

typedef struct mzid_pep mzid_pep;
struct mzid_pep {
	char id[100];
	char sequence[100];
	char mods[100];
};

typedef struct pep_info pep_info;
struct pep_info {
		char peptide[100];
		char mods[512];
		char **accessions;
		int accessions_cnt;
		char spectrum_id[100];
		double mz;
		double score;
		int charge;
		double rt;
		int *diags;
		double ion_mobility;
		int printed;
		};

typedef struct mgf_map mgf_map;
struct mgf_map {
		pep_info *peps;
		char mgf_name[100];
		char mgf_ref[100];
		int pep_cnt;
};

typedef struct msms_data msms_data;
struct msms_data {
                double mass;
                double intensity;
                int charge;
                };

typedef struct mgf_data mgf_data;
struct mgf_data {
                char title[512];
                double rt;
                double precursor_mass;
                int charge;
                msms_data *msms;
                int msms_values_cnt;
                };
typedef enum {
	SCIEX,
	BRUKER,
} mode_type;

void print_usage();
mgf_map *read_mzid(FILE *, int *, int);
char *move_ptr(char *,char);
char *skip_fields(char *, char, int);
void check_diag(pep_info **, int, double *, int, mgf_data *, int, mode_type); 
mgf_data *read_mgf_data(FILE *,int);
int count_spectra(FILE *);
int is_match(double, double, double);
void print_line(FILE *, char *,pep_info, int, mode_type);
/*penicillinFinder: looks for diagnostic ions in spectra of psms assigned by
PEAKS. Takes PEAKS search results for Bruker or SCIEX MS data in mzIdentML format plus mgf files.
By default, looks for 3 ions characteristic of penicillin haptenation, but other/selected 
diagnostic ions can be specified using -i option
*/
/*******************************************************************************/
void main(int argc,char **argv)
{
	FILE *f, *g;
	char mgf_fn[612], pep_fn[512], *pnt = NULL, *ptr = NULL, c, output[105], stem[100], filepath[512], diags[100], tmp[512];
	int mgf_cnt = 0, spectra_cnt = 0, i, j, k, m, diag_cnt = 0;
	mode_type mode = SCIEX;
	mgf_map *data = NULL;
	mgf_data *spectra = NULL;
	double *diagnostics, default_diag[3] = {160.0432,217.0647,335.1066};
	//char **proteins = NULL;

	pep_fn[0] = output[0] = stem[0] = diags[0] = '\0';
	strcpy(stem, "output");
	while ((c = getopt (argc, argv, "i:p:o:bh")) != -1) {
                switch (c) {
                        case 'b': //Mode BRUKER
				mode = BRUKER;
                                break;
                        case 'p':
                                strcpy(pep_fn, optarg);
                                break;
                        case 'o':
                                strcpy(stem, optarg);
                                break;
			case 'i':
				strcpy(diags, optarg);
				ptr = diags;
				while ((ptr = strchr(ptr, ',')) != NULL) {
					diag_cnt++;
					ptr += 1;
				}
				diag_cnt++;
				printf("%d diagnostic ion specified\n", diag_cnt);
        			if((diagnostics = calloc(diag_cnt,sizeof(double))) == NULL) {
                			printf("memory allocation error in main()\n");
                			exit(0);
        			}
				ptr = diags;
				for (i = 0; i < diag_cnt - 1; ++i) {
					pnt = move_ptr(ptr, ',');
					diagnostics[i] = atof(ptr);
					*pnt = ',';
					ptr = pnt + 1;
				}
				diagnostics[i] = atof(ptr);
				break;
			case 'h':
				print_usage();
				exit(0);
                        default: /*'?'*/			
				print_usage();
                                exit(0);
                }
        }
	if (!diag_cnt) {
		diag_cnt = 3;
        	if((diagnostics = calloc(diag_cnt,sizeof(double))) == NULL) {
                	printf("memory allocation error in main()\n");
                	exit(0);
        	}
		for (i = 0; i < diag_cnt; i++) {
			diagnostics[i] = default_diag[i];
		}
	}
	if (pep_fn[0] == '\0') {
		print_usage();
		exit(0);
	}
	printf("Running penicillinFinder on %s\n", pep_fn);
	for (i =0; i < diag_cnt; ++i)
		printf("ion %d: %lf\n", i, diagnostics[i]); 
//find filepath (Linux or Windows)
	filepath[0] = tmp[0] = '\0';
	strcpy(filepath,pep_fn);
	if (((pnt = strrchr(filepath, '\\')) != NULL) || ((pnt = strrchr(filepath, '/')) != NULL)) {
		pnt += 1;
		*pnt = '\0';
		strcpy(tmp,filepath);
		strcat(tmp,stem);
		strcpy(stem,tmp);
	}

	sprintf(output, "%s.csv", stem);
	if ((f = fopen(pep_fn,"r")) == NULL) {
               	printf("Can't open file %s\n",pep_fn);
               	exit(0);
        }
//store mzid data - psm info sorted by mgf
	data = read_mzid(f, &mgf_cnt, diag_cnt);
	fclose(f);
//read mgfs, check for diagnostic ions for each psm
	for (i=0; i<mgf_cnt; ++i) {
		tmp[0] = mgf_fn[0] = '\0';
		strcpy(mgf_fn,data[i].mgf_name);
		if (filepath[0] != '\0') {
			strcpy(tmp,filepath);
			strcat(tmp,mgf_fn);
			strcpy(mgf_fn,tmp);
		}
		if ((f = fopen(mgf_fn,"r")) == NULL) {
               		printf("Can't open file %s\n",mgf_fn);
               		exit(0);
        	}
		spectra_cnt = count_spectra(f);
		spectra = read_mgf_data(f, spectra_cnt);
		fclose(f);
		printf("%s:\nSpectra %d\n",data[i].mgf_name,spectra_cnt);
		check_diag(&data[i].peps, data[i].pep_cnt, diagnostics, diag_cnt, spectra, spectra_cnt, mode);
		//free spectra
		for (j = 0; j < spectra_cnt; ++j) {
			if (spectra[j].msms != NULL) 
				free(spectra[j].msms);
		}
		if (spectra != NULL)
			free(spectra);
		spectra = NULL;
	}
//print out results		
	if((g = fopen(output,"w")) == NULL){ 
                printf("Can't open file %s\n",output);
                exit(0);
        }
	fprintf(g, "Peptide,Modifications,\"-10lgP\",m/z,charge,RT(sec),");
	for (i = 0; i < diag_cnt; ++i) {
		fprintf(g, "%.2lf?,",diagnostics[i]);
	}
	fprintf(g,"total diagnostic ions seen,mgf,spectrum id,");
	if (mode == BRUKER) 
		fprintf(g, "1/k0,");
	fprintf(g, "protein accessions\n");

	for (i = 0; i < mgf_cnt; ++i) {  //print psms with BP-mod
		for (j = 0; j < data[i].pep_cnt; ++j) {
			if (strstr(data[i].peps[j].mods, "Benzylpenicillin") != NULL) {
				data[i].peps[j].printed = 1;
				print_line(g,data[i].mgf_name,data[i].peps[j],diag_cnt,mode);
			}
		}
	}
	for (i = 0; i < mgf_cnt; ++i) {  //print psms with other mods
		for (j = 0; j < data[i].pep_cnt; ++j) {
			if ((!data[i].peps[j].printed) && data[i].peps[j].mods[0] != '\0') {
				data[i].peps[j].printed = 1;
				print_line(g,data[i].mgf_name,data[i].peps[j],diag_cnt,mode);
			}
		}
	}
	fprintf(g, "\n");

	for (i = 0; i < mgf_cnt; ++i) {
		for (j = 0; j < data[i].pep_cnt; ++j) {
			if (!data[i].peps[j].printed) {
				print_line(g,data[i].mgf_name,data[i].peps[j],diag_cnt,mode);
			}
		}
	}
	fclose(g);
	//free data
	for (i = 0; i < mgf_cnt; ++i) {
		for (j = 0; j < data[i].pep_cnt; ++j) {
			if (data[i].peps[j].diags != NULL)
				free(data[i].peps[j].diags);
			for (k = 0; k < data[i].peps[j].accessions_cnt; ++k) {
				if (data[i].peps[j].accessions[k] != NULL)
					free(data[i].peps[j].accessions[k]);
			}
			if (data[i].peps[j].accessions != NULL)
				free(data[i].peps[j].accessions);
		}
		if (data[i].peps != NULL)
			free(data[i].peps);
	}
	if (data != NULL)
		free(data);
	if (diagnostics != NULL)
		free(diagnostics);
	printf("Finished\n");
	fflush(stdout);
	exit(0);
}
/*******************************************************************************/

/*******************************************************************************/
void print_line(FILE *f, char *mgf,pep_info pep, int diag_cnt, mode_type mode)
{
	int m, diags_total = 0;

	fprintf(f,"%s,%s,%lf,%lf,%d,%lf,",pep.peptide,pep.mods,pep.score,pep.mz,pep.charge,pep.rt);
	for (m = 0; m < diag_cnt; ++m) {
		fprintf(f,"%c,",(pep.diags[m]) ? 'Y':'\0');
		diags_total+=pep.diags[m];
	}
	fprintf(f,"%d,",diags_total);
	fprintf(f,"%s,%s",mgf,pep.spectrum_id);
	if (mode==BRUKER) 
		fprintf(f,",%lf",pep.ion_mobility);
	if (pep.accessions_cnt) {
		fprintf(f,",");
		for (m = 0; m < pep.accessions_cnt; ++m)
			fprintf(f,"%s;",pep.accessions[m]);
	}
	fprintf(f,"\n");

	return;
}
/*******************************************************************************/

/*******************************************************************************/
mgf_data *read_mgf_data(FILE *f,int spectra_cnt)
{
        mgf_data *tmp = NULL;
        char line[MAXSTR];
        int x = 0;
        int msms_cnt = 0;
        msms_data m[MAXIONS];

        if((tmp = (mgf_data *)calloc(spectra_cnt,sizeof(mgf_data))) == NULL) {
                printf("memory allocation error in read_mgf_data()\n");
                exit(0);
        }

        while(fgets(line,MAXSTR - 1, f) != NULL) {
                if(strncmp(line,"TITLE=",6) == 0) {
                        strcpy(tmp[x].title, line + 6);
                        while ((tmp[x].title[strlen(tmp[x].title)-1] == '\n') || (tmp[x].title[strlen(tmp[x].title)-1] == '\r')) {

                                tmp[x].title[strlen(tmp[x].title)-1] = '\0';
                        }
                }
                if(strncmp(line,"CHARGE=",7) == 0) sscanf(line + 7,"%d",&tmp[x].charge);
                if(strncmp(line,"PEPMASS=",8) == 0) sscanf(line + 8,"%lf",&tmp[x].precursor_mass);
                if(strncmp(line,"RTINSECONDS=",12) == 0) sscanf(line + 12,"%lf",&tmp[x].rt);
                if(isdigit(line[0])) {
                        sscanf(line,"%lf %lf %d",&m[msms_cnt].mass,&m[msms_cnt].intensity,&m[msms_cnt].charge);
                        msms_cnt ++;
                        if(msms_cnt >= MAXIONS) {
                                printf("Error - increase the size f MAXIONS\n");
                                exit(0);
                        }
                }
                if(strncmp(line,"END IONS",8) == 0) {
                        if((tmp[x].msms = (msms_data *)calloc(msms_cnt,sizeof(msms_data))) == NULL) {
                                printf("Memory allocation error\n");
                                exit(0);
                        }
                        memcpy(tmp[x].msms,m,sizeof(msms_data) * msms_cnt);
                        tmp[x].msms_values_cnt = msms_cnt;

                        msms_cnt = 0;

                        x ++;
                }
        }
 	return(tmp);
}
/*******************************************************************************/

/*******************************************************************************/
int count_spectra(FILE *f)
{
        char line[MAXSTR];
        int tmp = 0;

        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if(strncmp(line,"BEGIN IONS",10) == 0) {
                        tmp ++;
                }
        }
        rewind(f);
        return(tmp);
}
/*******************************************************************************/

/*******************************************************************************/
int is_match(double val1, double val2, double tol)
{
        if (isless((val1 - val2), tol) && isgreater((val1 - val2), -tol))
                return 1;
        else
                return 0;
}
/*******************************************************************************/

/*******************************************************************************/
void check_diag(pep_info **peps, int pep_cnt, double *diagnostics, int diag_cnt, mgf_data *spectra, int spectra_cnt, mode_type mode) 
{
	int i, j, k, m, match = 0;
	char match_title[105], *ptr = NULL, *pnt = NULL;

	for (i = 0; i < pep_cnt; ++i) {
		match_title[0] = '\0';
		if (mode == BRUKER) {
			sprintf(match_title, "%s,", (*peps)[i].spectrum_id);
		} 
		for (j = 0; j < spectra_cnt; ++j) {
			match = 0;
			if (mode == BRUKER) {
				if (strstr(spectra[j].title, match_title) != NULL) {
					match = 1;
				}
			} else {
				if (strcmp(spectra[j].title, (*peps)[i].spectrum_id) == 0) {
					match = 1;
				}
			}
			if (match) {
				if (i % 10000 == 0)
					printf("So far, matched %d peptides to spectra...\n",i); 
				(*peps)[i].rt = spectra[j].rt;
				if (mode == BRUKER) {
					if ((ptr = strstr(spectra[j].title, "IonMobility")) != NULL) {
						ptr += 13;
						pnt = move_ptr(ptr, '"');
						(*peps)[i].ion_mobility = atof(ptr);
						*pnt = '"';
					}
				}				
				for (k = 0; k < spectra[j].msms_values_cnt; ++k) {
					for (m = 0; m < diag_cnt; ++m) {
                				if (is_match(spectra[j].msms[k].mass, diagnostics[m], FRAGTOL))
                			       		(*peps)[i].diags[m] = 1;
					}
				}
			}
		}
	}
	return;
} 
/*******************************************************************************/

/*******************************************************************************/
void print_usage()
{
	printf("Usage:\nRequired\n");
	printf("\t-p\tSpecify peptide search results (peptides_1_1_0.mzid from PEAKS,\n");
	printf("\t\t\tin the same folder as the associated mgf files)\n");
	printf("Optional\n\t-b\tSpecify mode is Bruker to include ion mobility information\n");
	printf("\t-i\tSpecify diagnostic ions in comma-separated list\n\t\t(default is '160.0432,217.0647,335.1066')\n");
	printf("\t-o\tSpecify name for output csv file (default is 'output')\n");
	printf("\t-h\tPrint help\n");
	return;
}
/*******************************************************************************/

/*******************************************************************************/
mgf_map *read_mzid(FILE *f, int *mgf_cnt, int diag_cnt)
{
	char line[MAXSTR], *pnt = NULL, *ptr = NULL, mod[100], loc[5], modstring[100],res[20];
	char spectra_ref[100], title[100], pep_ref[100], **proteins, accession[100];
	int i,j, psm_cnt = 0, pep_cnt = 0, tmp_num = 0, pep_num = 0, index = -1, tmp_id = 0, maxprotid = 0, maxpepid = 0, pep_prot_index = 0;
	int *tmp_realloc = NULL;
	mzid_pep *peps;
	mgf_map *tmp;
	pep_info *tmp_map_peps;	
	peptide_evidence *pep_prot = NULL;
	
	modstring[0] = spectra_ref[0] = title[0] = '\0';
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (strstr(line, "<SpectraData ") != NULL) {
			*mgf_cnt +=1;
		}
		if (strstr(line, "<Peptide ") != NULL) {
			pep_cnt++;
			if ((ptr = strstr(line, "id=")) != NULL) {
				ptr += 12;
				pnt = move_ptr(ptr, '"');
				tmp_id = atoi(ptr);
				if (tmp_id > maxpepid)
					maxpepid = tmp_id;
				*pnt = '"';
			}
		}
		if (strstr(line, "<SpectrumIdentificationItem ") != NULL) {
			psm_cnt++;
		}
		if (strstr(line, "<DBSequence ") != NULL) {
			if ((ptr = strstr(line, " id=")) != NULL) {
				ptr += 16;
				pnt = move_ptr(ptr, '"');
				tmp_id = atoi(ptr);
				if (tmp_id > maxprotid)
					maxprotid = tmp_id;
				*pnt = '"';
			}
		}
	}
        rewind(f);
	printf("no. mgfs %d; peptides %d; psms %d\n", *mgf_cnt, pep_cnt, psm_cnt);
	fflush(stdout);
//get protein accessions
	if ((proteins = calloc(maxprotid + 1, sizeof(char *))) == NULL) {
               	printf("Memory allocation error in read_mzid()\n");
               	exit(0);
	}
	if ((pep_prot = calloc(maxpepid + 1, sizeof(peptide_evidence))) == NULL) {
               	printf("Memory allocation error in read_mzid()\n");
               	exit(0);
	}
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		//get accessions for each DBSequence
		accession[0] = '\0';
		if (strstr(line, "<DBSequence ") != NULL) {
			if ((ptr = strstr(line, "accession=")) != NULL) {
				ptr += 11;
				pnt = move_ptr(ptr, '"');
				strcpy(accession,ptr);
				*pnt = '"';
			}
			if (accession[0] == '\0') {
				printf("Couldn't find accession in DBSequence line %s\n",line);
				exit(0);
			}
			if ((ptr = strstr(line, "id=")) != NULL) {
				ptr += 15;
				pnt = move_ptr(ptr, '"');
				tmp_id = atoi(ptr);
				*pnt = '"';
				
				if ((proteins[tmp_id] = calloc(strlen(accession) + 1, sizeof(char))) == NULL) {
        			       	printf("Memory allocation error in read_mzid()\n");
        			       	exit(0);
				}
				strcpy(proteins[tmp_id],accession);
			}
		}
		//get all dbseq refs (i.e. index) for each peptide
		if (strstr(line, "<PeptideEvidence ") != NULL) {
			if ((ptr = strstr(line, "peptide_ref=")) != NULL) {
				pep_prot_index = 0;
				ptr += 21;
				pnt = move_ptr(ptr, '"');
				pep_prot_index = atoi(ptr);
				*pnt = '"';
			}
			if (!pep_prot_index) {
				printf("Couldn't find peptide id in Peptide Evidence line %s\n",line);
				exit(0);
			}
			if ((ptr = strstr(line, "dBSequence_ref=")) != NULL) {
				ptr += 27;
				pnt = move_ptr(ptr, '"');
				tmp_id = atoi(ptr);
				*pnt = '"';
				if (!pep_prot[pep_prot_index].cnt) {	
					if ((pep_prot[pep_prot_index].prot_ids = calloc(1, sizeof(int))) == NULL) {
        				       	printf("Memory allocation error in read_mzid()\n");
        			       		exit(0);
					}
				} else {
					if ((tmp_realloc = realloc(pep_prot[pep_prot_index].prot_ids, ((pep_prot[pep_prot_index].cnt + 1) * sizeof(int)))) == NULL) {
        				       	printf("Memory allocation error in read_mzid()\n");
                        			exit(0);
                			} else {
                        			pep_prot[pep_prot_index].prot_ids = tmp_realloc;
					}
				}
				pep_prot[pep_prot_index].prot_ids[pep_prot[pep_prot_index].cnt] = tmp_id;
				pep_prot[pep_prot_index].cnt++;
			}
		}
	}
	rewind(f);

	if ((peps = calloc(pep_cnt, sizeof(mzid_pep))) == NULL) {
               	printf("Memory allocation error in read_mzid()\n");
               	exit(0);
	}
	if ((tmp = calloc(*mgf_cnt, sizeof(mgf_map))) == NULL) {
               	printf("Memory allocation error in read_mzid()\n");
               	exit(0);
	}
	for (i=0; i < *mgf_cnt; ++i) {
		if ((tmp[i].peps = calloc(psm_cnt, sizeof(pep_info))) == NULL) {
               		printf("Memory allocation error in read_mzid()\n");
               		exit(0);
		}
	}

	while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (strstr(line, "<SpectraData ") != NULL) {
			if ((ptr = strstr(line, "location=")) != NULL) {
				ptr += 10;
				while ((pnt = strchr(ptr, '/')) != NULL) {
					ptr = pnt +1;
				}
				pnt = move_ptr(ptr, '"');
				strcpy(tmp[tmp_num].mgf_name, ptr);
				*pnt = '"';
			}
			if ((ptr = strstr(line, "id=")) != NULL) {
				ptr += 4;
				pnt = move_ptr(ptr,'"');
				strcpy(tmp[tmp_num].mgf_ref, ptr);
				*pnt = '"';
			}
			tmp_num++;
		}
		if (strstr(line, "<Peptide ") != NULL) {
			if ((ptr = strstr(line, "id=")) != NULL) {
				ptr += 4;
				pnt = move_ptr(ptr, '"');
				strcpy(peps[pep_num].id ,ptr);
				*pnt = '"';				
			}
			while(fgets(line, MAXSTR-1, f) != NULL && strstr(line,"</Peptide>") == NULL) {
				if ((ptr = strstr(line, "<PeptideSequence>")) != NULL) {
					ptr += 17;
					pnt = move_ptr(ptr, '<');
					strcpy(peps[pep_num].sequence ,ptr);
					*pnt = '<';				
				}
				if (strstr(line, "<Modification ") != NULL) {
					mod[0] = loc[0] = res[0] = '\0';
					
					if ((ptr = strstr(line, "location=")) != NULL) {
						ptr += 10;
						pnt = move_ptr(ptr, '"');
						strcpy(loc ,ptr);
						*pnt = '"';				
					}
					if ((ptr = strstr(line, "residues=")) != NULL) {
						ptr += 10;
						pnt = move_ptr(ptr, '"');
						strcpy(res ,ptr);
						*pnt = '"';				
					}
					while(fgets(line, MAXSTR-1, f) != NULL && strstr(line,"</Modification>") == NULL) {
						if (strstr(line, "<cvParam ") != NULL) {
							if ((ptr = strstr(line, "name=")) != NULL) {
								ptr += 6;
								pnt = move_ptr(ptr, '"');
								strcpy(mod,ptr);
								*pnt = '"';
								
								if (strcmp(mod,"unknown modification") == 0) {
									mod[0] = '\0';
									if ((ptr = strstr(line, "value=")) != NULL) {
										ptr += 7;
										pnt = move_ptr(ptr, '"');
										strcpy(mod,ptr);
										*pnt = '"';
									}	
								}			
							}
						}
					}
					sprintf(modstring, "%s%s(%s)@%s;",modstring,mod,res,loc);
				}
			}
			strcpy(peps[pep_num].mods, modstring);
			modstring[0] = '\0';
			pep_num++;
		}
	}
        rewind(f);
	while(fgets(line,MAXSTR - 1,f) != NULL) {
		if (strstr(line, "<SpectrumIdentificationResult ") != NULL) {
			index = -1;
			if ((ptr = strstr(line, "spectraData_ref=")) != NULL) {
				ptr += 17;
				pnt = move_ptr(ptr, '"');
				strcpy(spectra_ref, ptr);
				*pnt = '"';
			}
			for (i = 0; i < *mgf_cnt; ++i) {
				if (strcmp(tmp[i].mgf_ref, spectra_ref) == 0) {
					index = i;
				}
			}
			if (index < 0) {
				printf("Couldn't match spectra_ref\n");
				exit(0);
			}
			if ((ptr = strstr(line, "spectrumID=")) != NULL) {
				ptr += 12;
				pnt = move_ptr(ptr,'"');
				strcpy(title, ptr);
				*pnt = '"';
			}
			while(fgets(line,MAXSTR-1,f) != NULL && strstr(line,"</SpectrumIdentificationResult>") ==NULL) {
				if (strstr(line,"<SpectrumIdentificationItem ") != NULL) {
					pep_ref[0] = '\0';
					strcpy(tmp[index].peps[tmp[index].pep_cnt].spectrum_id, title);
					if ((ptr = strstr(line,"chargeState=")) != NULL) {
						ptr += 13;
						pnt = move_ptr(ptr, '"');
						tmp[index].peps[tmp[index].pep_cnt].charge = atoi(ptr);
						*pnt = '"';
					}
					if ((ptr = strstr(line,"experimentalMassToCharge=")) != NULL) {
						ptr += 26;
						pnt = move_ptr(ptr, '"');
						tmp[index].peps[tmp[index].pep_cnt].mz = atof(ptr);
						*pnt = '"';
					}
					if ((ptr = strstr(line,"peptide_ref=")) != NULL) {
						ptr += 13;
						pnt = move_ptr(ptr, '"');
						strcpy(pep_ref,ptr);
						*pnt = '"';
						
						for (i=0; i<pep_cnt; ++i) {
							if (strcmp(peps[i].id,pep_ref) == 0) {
								strcpy(tmp[index].peps[tmp[index].pep_cnt].peptide, peps[i].sequence);
								strcpy(tmp[index].peps[tmp[index].pep_cnt].mods, peps[i].mods);
								ptr = pep_ref;
								ptr += 8;
								pep_prot_index = atoi(ptr);
								if (!pep_prot_index) {
									printf("not finding index in %s\n", peps[i].id);
									exit(0);
								}
								//get accession refs, copy accessions to tmp[index].peps[tmp[index].pep_cnt].accessions
								if (pep_prot[pep_prot_index].cnt) {
									tmp[index].peps[tmp[index].pep_cnt].accessions_cnt = pep_prot[pep_prot_index].cnt;
									if((tmp[index].peps[tmp[index].pep_cnt].accessions = calloc(pep_prot[pep_prot_index].cnt,sizeof(char *))) == NULL) {
                								printf("Memory allocation error in read_mzid()\n");
               								        exit(0);
                							}
									for (i = 0; i < pep_prot[pep_prot_index].cnt; ++i) {
										if((tmp[index].peps[tmp[index].pep_cnt].accessions[i] = calloc(strlen(proteins[pep_prot[pep_prot_index].prot_ids[i]]) + 1,sizeof(char))) == NULL) {
                									printf("Memory allocation error in read_mzid()\n");
               								        	exit(0);
                								}
										strcpy(tmp[index].peps[tmp[index].pep_cnt].accessions[i],proteins[pep_prot[pep_prot_index].prot_ids[i]]);
									}
								}
								break;
							}
						}
					}
					while(fgets(line,MAXSTR-1,f) != NULL && strstr(line,"</SpectrumIdentificationItem>") ==NULL) {
						if (strstr(line,"<cvParam ") != NULL) {
							if ((ptr = strstr(line,"value=")) != NULL) {
								ptr += 7;
								pnt = move_ptr(ptr, '"');
								tmp[index].peps[tmp[index].pep_cnt].score = atof(ptr);
								*pnt = '"';
							}
						}
					}
					tmp[index].pep_cnt++;
				}
			}
		}
	} 
	//copy to peps to correct sized array and free other
        for (i=0; i < *mgf_cnt; ++i) {
		if((tmp_map_peps = calloc(tmp[i].pep_cnt,sizeof(pep_info))) == NULL) {
                	printf("Memory allocation error in read_mzid()\n");
                        exit(0);
                }
                memcpy(tmp_map_peps,tmp[i].peps,sizeof(pep_info) * tmp[i].pep_cnt);
		if (tmp[i].peps != NULL)
			free(tmp[i].peps);
		tmp[i].peps = tmp_map_peps;
		tmp_map_peps = NULL;
	}
	for (i = 0; i < *mgf_cnt; ++i) {
		for (j = 0; j < tmp[i].pep_cnt; ++j) {
			if ((tmp[i].peps[j].diags = calloc(diag_cnt, sizeof(int))) == NULL) {
        	       		printf("Memory allocation error in read_mzid()\n");
        	       		exit(0);
			}
		}
	}
//clean up
	for (i = 0; i < maxprotid + 1; ++i) {
		if (proteins[i] != NULL)
			free(proteins[i]);
	}
	if (proteins != NULL)
		free(proteins);
	for (i = 0; i < maxpepid + 1; ++i) {
		if (pep_prot[i].prot_ids != NULL)
			free(pep_prot[i].prot_ids);
	}
	if (pep_prot != NULL)
		free(pep_prot);
	if (peps != NULL)
		free(peps);

	return tmp;
}
/*******************************************************************************/

/*******************************************************************************/
char *move_ptr(char *str,char delim)
{
	char *tmp = NULL;

	if((tmp = strchr(str,delim)) != NULL) *tmp = '\0';

	return(tmp);
}
/*******************************************************************************/

/*******************************************************************************/
char *skip_fields(char *str, char delim, int num)
{
	char *tmp = NULL, *ptr;
	int i;
	
	ptr = str;
	for (i = 0; i < num; ++i) {
		if ((tmp = strchr(ptr,delim)) != NULL && *(tmp + 1) != '\0') {
			ptr = tmp +1;
		} else {
			printf("Can't skip %d fields\n", num);
		}
	}

	return(ptr);
}
/*******************************************************************************/

/*******************************************************************************/


