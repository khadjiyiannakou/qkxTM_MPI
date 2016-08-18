#include <string.h> 
#include <math.h>
#include <stdarg.h>
extern "C"{
#include <lime.h>
}
#include <lattice_util.h>

using namespace qkxTM;

static int qkxTM_isBigEndian()
{
  union{
    char C[4];
    int  R   ;
  }word;
  word.R=1;
  if(word.C[3]==1) return 1;
  if(word.C[0]==1) return 0;
  return -1;
}


static void qkxTM_swap_8(double *Rd, long int N)
{
  register char *i,*j,*k;
  char swap;
  char *max;
  char *R = (char*) Rd;

  max = R+(N<<3);
  for(i=R;i<max;i+=8)
    {
      j=i; k=j+7;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
    }
}

/*
static char* qkxTM_getParams(char* fname,int *len)
{
  FILE *pfile;
  char *params;
  int i;
  size_t size;

  if ((pfile=fopen(fname,"r"))==NULL)
    {
      fprintf(stderr,"Error, cannot open %s for reading\n",fname);
      return(NULL);
    }

  i=0;
  while(!feof(pfile))
    {
      fgetc(pfile);
      i++;
    }
  *(len)=i;
  rewind(pfile);
  params = (char*)malloc(*(len)*sizeof(char));

  size = fread(params,sizeof(char),*len,pfile);

  fclose(pfile);

  return params;
}
*/

static char* qkxTM_getParam(char token[],char* params,int len)
{
  int i,token_len=strlen(token);

  for(i=0;i<len-token_len;i++)
    {
      if(memcmp(token,params+i,token_len)==0)
	{
	  i+=token_len;
	  *(strchr(params+i,'<'))='\0';
	  break;
	}
    }
  return params+i;
}

/*
static char* qkxTM_getParamComma(char token[],char* params,int len)
{
  int i,token_len=strlen(token);

  for(i=0;i<len-token_len;i++)
    {
      if(memcmp(token,params+i,token_len)==0)
        {
          i+=token_len;
          *(strchr(params+i,','))='\0';
          break;
        }
    }
  return params+i;
}
*/

int qkxTM_MPI_getGaugeLime(char *fname, Gauge &u)
{

   LimeReader *limereader;
   FILE *fid;
   char *lime_type,*lime_data;
   n_uint64_t lime_data_size;
   char dummy;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   unsigned long int i=0;//,j=0;
   //   unsigned int stride;
   unsigned short chunksize; //,mu,nu,c1,c2;
   //   char *swap[16*3*3];
   char *buffer;
   //   unsigned int x,y,z,t;
   //   int  isDouble;
   // int error_occured=0;
   int rank = u.rank;
   int precision;
   bool isDouble;
   
     
   if(u.initialized == false) errorQKXTM("gauge not initialized\n");
   
   if(rank == 0)
   {
      /* read lime header */
      fid=fopen(fname,"r");
      if(fid==NULL) errorQKXTM("cant open configuration for reading\n");
      if ((limereader = limeCreateReader(fid))==NULL) errorQKXTM("Error cant create limeReader\n");

      while(limeReaderNextRecord(limereader) != LIME_EOF )
	{
	  lime_type = limeReaderType(limereader);
	  if(strcmp(lime_type,"ildg-binary-data")==0)
            {
	      break;
            }
	  if(strcmp(lime_type,"ildg-format")==0)
            {
	      lime_data_size = limeReaderBytes(limereader);
	      lime_data = (char * )malloc(lime_data_size);
	      limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
	      sscanf(qkxTM_getParam((char*)"<precision>",lime_data, lime_data_size),"%i",&precision);    
	      //printf("got precision: %i\n",isDouble);
	      free(lime_data);
            }
	}
      /* read 1 byte to set file-pointer to start of binary data */
      lime_data_size=1;
      limeReaderReadData(&dummy,&lime_data_size,limereader);
      offset = ftell(fid)-1;
      limeDestroyReader(limereader);      
      fclose(fid);
      
   }// if root condition
      
   MPI_Bcast(&precision,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);

   if(precision == 64)
      isDouble=true;      
   else
   {
     errorQKXTM("Unknown precision in configuration\n");
   }   
   

   if(isDouble == true)
   {
      sizes[0]=u.L[3]; 
      sizes[1]=u.L[2];
      sizes[2]=u.L[1];
      sizes[3]=u.L[0];
      sizes[4]=NLINKS;
      lsizes[0]=u.lL[3];
      lsizes[1]=u.lL[2];
      lsizes[2]=u.lL[1];
      lsizes[3]=u.lL[0];
      lsizes[4]=sizes[4];
      starts[0]=u.procPos[3]*lsizes[0];
      starts[1]=u.procPos[2]*lsizes[1];
      starts[2]=u.procPos[1]*lsizes[2];
      starts[3]=u.procPos[0]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, (char*)"native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:
      chunksize=NLINKS*sizeof(double);
      size_t memoryN = chunksize * ((size_t)u.lV4d);
      buffer = (char*) malloc(memoryN);

   //      printf("%d\n",chunksize*u.lV4d);
   //      fflush(stdout);


      if(buffer==NULL) errorQKXTM("Error cant allocate memory for buffer to read Configuration\n");



      int error_code = MPI_File_read_all(mpifid, buffer, NLINKS * u.lV4d, MPI_DOUBLE, &status);
      if(error_code != MPI_SUCCESS){
	printf("problem\n");
	fflush(stdout);
      }


      if(!qkxTM_isBigEndian()){
	qkxTM_swap_8((double*) buffer,NLINKS * u.lV4d);
      }

      i=0;
      

      for(int t=0; t < u.lL[3]; t++)
	for(int z=0; z < u.lL[2]; z++)
	  for(int y=0; y < u.lL[1]; y++)
	    for(int x=0; x < u.lL[0]; x++)
	      for(int mu=0; mu < 4; mu++)
		{
		  memcpy(&(u.M[LEXIC(t,z,y,x,u.lL)][mu][0][0].real()),&(buffer[i]),144);
		  i+=144;
		}

      // for(int t=0; t < u.lL[3]; t++)
      // 	for(int z=0; z < u.lL[2]; z++)
      // 	  for(int y=0; y < u.lL[1]; y++)
      // 	    for(int x=0; x < u.lL[0]; x++)
      // 	      for(int mu=0; mu < 4; mu++)
      // 		for(int c1 = 0 ; c1 < 3 ; c1++)
      // 		  for(int c2 = 0 ; c2 < 3 ; c2++)
      // 		    {
      // 		      if(isnan(u.M[LEXIC(t,z,y,x,u.lL)][mu][c1][c2].real())){
      // 			printf("nan in %d,%d,%d,%d \n",t,z,y,x);
      // 			fflush(stdout);
      // 			exit(-1);
      // 		      }
      // 		    }

      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);
      
      return 0;
   }//end isDouble condition

   /*
   else
   {
      sizes[0]=u.L[0]; 
      sizes[1]=u.L[3];
      sizes[2]=u.L[2];
      sizes[3]=u.L[1];
      sizes[4]=NLINKS;
      lsizes[0]=u.lL[0];
      lsizes[1]=u.lL[3];
      lsizes[2]=u.lL[2];
      lsizes[3]=u.lL[1];
      lsizes[4]=sizes[4];
      starts[0]=u.procPos[0]*lsizes[0];
      starts[1]=u.procPos[3]*lsizes[1];
      starts[2]=u.procPos[2]*lsizes[2];
      starts[3]=u.procPos[1]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:

      chunksize=NLINKS*sizeof(float);

      buffer = (char*) malloc(chunksize*u.lV4d);
      if(buffer==NULL)
      {
         fprintf(stderr,"process %i: Error in qkxTM_getGaugeLime! Out of memory\n",u.myid);
         return(1);
      }
      MPI_File_read_all(mpifid, buffer, 4*3*3*2*u.lV, MPI_FLOAT, &status);
     
      if(!qkxTM_isBigEndian())      
         qkxTM_swap_4((float*) buffer,2*4*3*3*u.lV3);
      
      for(t=0; t<u.lL[0];t++)
      for(z=0; z<u.lL[3];z++)
      for(y=0; y<u.lL[2];y++)
      for(x=0; x<u.lL[1];x++)
      for(mu=0; mu<4; mu++)
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
      {
         nu=(mu+1)%4;
         j=qkxTM_LEXIC(t,x,y,z,u.lL);
         u->D[j][nu][c1][c2].re = *((float*) &(buffer[i])); i+=4;
         u->D[j][nu][c1][c2].im = *((float*) &(buffer[i])); i+=4;
      }      
            
      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);            

      return 0;
   }//end isDouble condition

   */
   return 0;
}
 
