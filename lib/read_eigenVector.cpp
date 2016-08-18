#include <string.h> 
#include <math.h>
#include <stdarg.h>
extern "C"{
#include <lime.h>
}
#include <lattice_util.h>

using namespace qkxTM;

static void qcd_swap_4(float *Rd, size_t N)
{
  register char *i,*j,*k;
  char swap;
  char *max;
  char *R =(char*) Rd;

  max = R+(N<<2);
  for(i=R;i<max;i+=4)
    {
      j=i; k=j+3;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
    }
}


static void qcd_swap_8(double *Rd, int N)
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

static int qcd_isBigEndian()
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

static char* qcd_getParam(char token[],char* params,int len)
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

void qkxTM_MPI_readEigenVectors(char *fname, Vector &v){

   LimeReader *limereader;
   FILE *fid;
   char *lime_type,*lime_data;
   unsigned long int lime_data_size;
   char dummy;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   unsigned int i,j;
   unsigned short int chunksize,mu,c1;
   char *buffer;
   unsigned int x,y,z,t;
   int  isDouble;
   int error_occured=0;
   int next_rec_is_prop = 0;
   
   v.zero();

   if(v.rank == 0)
     {
       /* read lime header */
       fid=fopen(fname,"r");
       if(fid==NULL)
	 {
	   fprintf(stderr,"process 0: Error in %s Could not open %s for reading\n",__func__, fname);
	   error_occured=1;
	 }
       if ((limereader = limeCreateReader(fid))==NULL)
	 {
	   fprintf(stderr,"process 0: Error in %s! Could not create limeReader\n", __func__);
	   error_occured=1;
	 }
       if(!error_occured)
	 {
	   while(limeReaderNextRecord(limereader) != LIME_EOF )
	     {
	       lime_type = limeReaderType(limereader);
	       if(strcmp(lime_type,"propagator-type")==0)
		 {
		   lime_data_size = limeReaderBytes(limereader);
		   lime_data = (char * )malloc(lime_data_size);
		   limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
		   
		   if (strncmp ("DiracFermion_Source_Sink_Pairs", lime_data, strlen ("DiracFermion_Source_Sink_Pairs"))!=0 &&
		       strncmp ("DiracFermion_Sink", lime_data, strlen ("DiracFermion_Sink"))!=0 )
		     {
		       fprintf (stderr, " process 0: Error in %s! Got %s for \"propagator-type\", expecting %s or %s\n", __func__, lime_data, 
				"DiracFermion_Source_Sink_Pairs", 
				"DiracFermion_Sink");
		       error_occured = 1;
		       break;
		     }
		   free(lime_data);
		 }
	       //lime_type="scidac-binary-data";
	       if((strcmp(lime_type,"etmc-propagator-format")==0) || (strcmp(lime_type,"etmc-source-format")==0) || (strcmp(lime_type,"etmc-eigenvectors-format")==0) )
		 {
		   lime_data_size = limeReaderBytes(limereader);
		   lime_data = (char * )malloc(lime_data_size);
		   limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
		   sscanf(qcd_getParam("<precision>",lime_data, lime_data_size),"%i",&isDouble);    
		   //printf("got precision: %i\n",isDouble);
		   free(lime_data);
		   
		   next_rec_is_prop = 1;
		 }
	       if(strcmp(lime_type,"scidac-binary-data")==0 && next_rec_is_prop)
		 {	      
		   break;
		 }
	     }
	   /* read 1 byte to set file-pointer to start of binary data */
	   lime_data_size=1;
	   limeReaderReadData(&dummy,&lime_data_size,limereader);
	   offset = ftell(fid)-1;
	   limeDestroyReader(limereader);      
	   fclose(fid);
	 }     
     }//end myid==0 condition
   
   MPI_Bcast(&error_occured,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_occured) return;
   
   MPI_Bcast(&isDouble,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);


   
   if(isDouble==64)
     isDouble=1;      
   else if(isDouble==32)
     isDouble=0; 
   else
     {
       fprintf(stderr,"process %i: Error in %s! Unsupported precision\n",v.rank, __func__);
     }  
   
   if(isDouble)
     {
       
       sizes[0]=v.L[3];
       sizes[1]=v.L[2];
       sizes[2]=v.L[1];
       sizes[3]=v.L[0];
       sizes[4]=4*3*2;
       lsizes[0]=v.lL[3];
       lsizes[1]=v.lL[2];
       lsizes[2]=v.lL[1];
       lsizes[3]=v.lL[0];
       lsizes[4]=sizes[4];
       starts[0]=v.procPos[3]*lsizes[0];
       starts[1]=v.procPos[2]*lsizes[1];
       starts[2]=v.procPos[1]*lsizes[2];
       starts[3]=v.procPos[0]*lsizes[3];
       starts[4]=0;
       
       MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
       MPI_Type_commit(&subblock);
       
       MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
       MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
       
       //load time-slice by time-slice:
       chunksize=4*3*2*sizeof(double);
       buffer = (char*) malloc(chunksize*v.lV4d);
       if(buffer==NULL)
	 {
	   fprintf(stderr,"process %i: Error in %s! Out of memory\n",v.rank, __func__);
	   return;
	 }
       MPI_File_read_all(mpifid, buffer, 4*3*2*v.lV4d, MPI_DOUBLE, &status);
       if(!qcd_isBigEndian())      
	 qcd_swap_8((double*) buffer,(size_t)(2*4*3)*(size_t)v.lV4d);
       i=0;



       for(t=0; t<v.lL[3];t++)
	 for(z=0; z<v.lL[2];z++)
	   for(y=0; y<v.lL[1];y++)
	     for(x=0; x<v.lL[0];x++)
	       for(mu=0; mu<4; mu++)
		 for(c1=0; c1<3; c1++)
		   {
		     
		     int oddBit     = (x+y+z+t) & 1;
		     if(oddBit){
		       v.M[LEXIC(t,z,y,x,v.lL)][mu][c1].real() = ((double*)buffer)[i];
		       v.M[LEXIC(t,z,y,x,v.lL)][mu][c1].imag() = ((double*)buffer)[i+1];
		       i+=2;
		     }
		     else{
		       i+=2;
		     }
		     
		   }

	   


       free(buffer);
       MPI_File_close(&mpifid);
       MPI_Type_free(&subblock);
       
       
     }//end isDouble condition
   else
     {
       sizes[0]=v.L[3];
       sizes[1]=v.L[2];
       sizes[2]=v.L[1];
       sizes[3]=v.L[0];
       sizes[4]=4*3*2;
       lsizes[0]=v.lL[3];
       lsizes[1]=v.lL[2];
       lsizes[2]=v.lL[1];
       lsizes[3]=v.lL[0];
       lsizes[4]=sizes[4];
       starts[0]=v.procPos[3]*lsizes[0];
       starts[1]=v.procPos[2]*lsizes[1];
       starts[2]=v.procPos[1]*lsizes[2];
       starts[3]=v.procPos[0]*lsizes[3];
       starts[4]=0;
       
       
       MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
       MPI_Type_commit(&subblock);
       
       MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
       MPI_File_set_view(mpifid, offset, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
       
       //load time-slice by time-slice:
       chunksize=4*3*2*sizeof(float);
       buffer = (char*) malloc(chunksize*v.lV4d);
       if(buffer==NULL)
	 {
	   fprintf(stderr,"process %i: Error in %s! Out of memory\n",v.rank, __func__);
	   return;
	 }
       MPI_File_read_all(mpifid, buffer, 4*3*2*v.lV4d, MPI_FLOAT, &status);
       
       if(!qcd_isBigEndian())
	 qcd_swap_4((float*) buffer,(size_t)(2*4*3)*(size_t)v.lV4d);
       
       i=0;


       for(t=0; t<v.lL[3];t++)
	 for(z=0; z<v.lL[2];z++)
	   for(y=0; y<v.lL[1];y++)
	     for(x=0; x<v.lL[0];x++)
	       for(mu=0; mu<4; mu++)
		 for(c1=0; c1<3; c1++)
		   {
		     int oddBit     = (x+y+z+t) & 1;
		     if(oddBit){
		       v.M[LEXIC(t,z,y,x,v.lL)][mu][c1].real() = *((float*)(buffer + i));
		       v.M[LEXIC(t,z,y,x,v.lL)][mu][c1].imag() = *((float*)(buffer + i + 4));
		       i+=8;
		     }
		     else{
		       i+=8;
		     }
		   }      

       
       
       free(buffer);
       MPI_File_close(&mpifid);
       MPI_Type_free(&subblock);            
       
       
     }//end isDouble condition
   
   
}//end qcd_getVectorLime 
