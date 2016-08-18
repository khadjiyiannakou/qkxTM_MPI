#include <string.h> 
#include <math.h>
#include <stdarg.h>
extern "C" {
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


static void qkxTM_swap_8(double *Rd, int N)
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




/*writes a double precision ildg formatted gauge configuration file*/
int qkxTM_writeGaugeLime(char *fname, Gauge &u, char* message)
{
   FILE *fid;   
   int error_in_header=0;
   LimeWriter *limewriter;
   LimeRecordHeader *limeheader;
   int ME_flag=0, MB_flag=0, limeStatus;
   unsigned long int message_length;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   unsigned long int i;
   int chunksize,mu,nu;
   char *buffer;
   int x,y,z,t;
   char ildgHeader[2048];
   int rank = u.rank;

   //process 0 creates file and writes header   
   if(rank == 0)
   {
      fid=fopen(fname,"w");
      if(fid==NULL)
      {
         fprintf(stderr,"process 0: Error in qkxTM_writeGaugeLime! Could not open %s for writing\n",fname);
         error_in_header=1;
      }else
      {

         limewriter = limeCreateWriter(fid); 
         if(limewriter == (LimeWriter*)NULL) {
            fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME error in file %s for writing!\n", fname);
            error_in_header=1;
         }else
         {            
            sprintf(ildgHeader, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqkxTM.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqkxTM.org/ildg filefmt.xsd\">\n  <version> 1.0 </version>\n  <field> su3gauge </field>\n  <precision> 64 </precision>\n  <lx> %d </lx>\n  <ly> %d </ly>\n  <lz> %d </lz>\n  <lt> %d </lt>\n</ildgFormat>",u.L[0],u.L[1], u.L[2], u.L[3]);
            message_length=(unsigned long int) strlen(ildgHeader);
            MB_flag=1; ME_flag=0;
            limeheader = limeCreateHeader(MB_flag, ME_flag,(char*) "ildg-format", message_length);
            if(limeheader == (LimeRecordHeader*)NULL) 
            {
               fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME create header ildg-format error.\n");
               error_in_header=1;
            }
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME write header ildg-format error %d\n", limeStatus);
              error_in_header=1;
            }
            limeDestroyHeader(limeheader);
            limeStatus = limeWriteRecordData(ildgHeader, &message_length, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME write header ildg-format error %d\n", limeStatus);
              error_in_header=1;
            }

            message_length=strlen(message);
            ME_flag=0; MB_flag=0;
            limeheader = limeCreateHeader(MB_flag, ME_flag,(char*) "xlf-info", message_length);
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
               fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME write header error %d\n", limeStatus);
               error_in_header=1;
            }
            limeDestroyHeader( limeheader );
            limeWriteRecordData(message, &message_length, limewriter);


            message_length = u.V4d*4*3*3*sizeof(Complex);
            MB_flag=0; ME_flag=1;
            limeheader = limeCreateHeader(MB_flag, ME_flag,(char*) "ildg-binary-data", message_length);
            limeStatus = limeWriteRecordHeader( limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
               fprintf(stderr, "Error in qkxTM_writeGaugeLime. LIME write header error %d\n", limeStatus);
               error_in_header=1;
            }
            limeDestroyHeader( limeheader );
         }
         //make one fake record-write to set offset
         message_length=1;
         limeWriteRecordData(message, &message_length, limewriter);
         offset = ftell(fid)-1;
         fclose(fid);
      }//fid non null condition
   }//myid==0 condition
   
   MPI_Bcast(&error_in_header,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_in_header)
      return(1);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);
   
   //header written. Now use MPI-I/O to write binary data
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

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, (char*)"native", MPI_INFO_NULL);

   chunksize=4*3*3*sizeof(Complex);
   buffer = (char*) malloc(chunksize*u.lV4d);
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qkxTM_writeGaugeLime! Out of memory\n",rank);
      return(1);
   }
   i=0;

   for(t=0; t < u.lL[3]; t++)
     for(z=0; z < u.lL[2]; z++)
       for(y=0; y < u.lL[1]; y++)
	 for(x=0; x < u.lL[0]; x++)
	   for(mu=0; mu < 4; mu++)
	     
	     {
	       memcpy(&(buffer[i]),&(u.M[LEXIC(t,z,y,x,u.lL)][mu][0][0].real()),144);
	       i+=144;
	     }

   if(!qkxTM_isBigEndian())
     qkxTM_swap_8((double*) buffer,(size_t)(2*4*3*3)*(size_t)u.lV4d);
   
   MPI_File_write_all(mpifid, buffer, 4*3*3*2*u.lV4d, MPI_DOUBLE, &status);

   free(buffer);
   MPI_File_close(&mpifid);
   MPI_Type_free(&subblock);

   return 0;
}//end qkxTM_writeGaugeLime
