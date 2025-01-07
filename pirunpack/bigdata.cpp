
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
//#include <sys/time.h>

#include "generic_utils.h"
#include "spqlios/lagrangehalfc_impl.h"
#include <chrono>
#include "64header.h"

#include "poc_64types.h"

Random* global_random = new Random();
const int Globals::k=1;
const int Globals::N=2048; 
const int Globals::t=12;
const int Globals::smalln=N; 
const int Globals::bgbit=15;
const int Globals::l=2;
const int Globals::basebit=4;
using namespace std;


Globals::Globals(){




        
        // Offset
        torusDecompOffset = 0;
        for (int i = 0; i <= l; ++i) torusDecompOffset |= (UINT64_C(1)<<(63-i*bgbit));
        
        // Buffer
        torusDecompBuf = new uint64_t[N]; 
        
        //secret keys (and their polynomial interpretation)

       lwekey = new int[N];
        for (int i = 0; i < N; ++i) lwekey[i] = random_bit();
        
        tlwekey = new IntPolynomiala(N);
        for (int i = 0; i < N; ++i) tlwekey->coefs[i] = lwekey[i]; 
        // level 2
        in_key = new int64_t[N+1];
        for (int i = 0; i < N; ++i) in_key[i] = lwekey[i];
        in_key[N] = -1; 
	
	//temp = new TLweSample64(N);
     

  int base = 1 <<basebit;

 privKS = new_array4<TLweSample64>(k+1, N+1, t, base, N);
        for (int z = 0; z <= k; ++z) {
            for (int i = 0; i <= N; ++i) {
                for (int j = 0; j < t; ++j) {
                    for (int u = 0; u < base; ++u) {
                        Torus64 messPrivKS = (in_key[i] << (64-(j+1)*basebit)) * u;
                        tLwe64EncryptZero(&privKS[z][i][j][u],pow(2.,-55), this);
			
                        privKS[z][i][j][u].a[z].coefs[0] += messPrivKS;
                    }
                }
            }
	}


  }

  int64_t array_max(Torus64Polynomial *a, int size)
  {
	

	int64_t max = abs(a->coefs[0]);

	for(int i=1;i< size;i++) {

	int64_t temp= abs(a->coefs[i]);
      if(temp > max) max=temp;
	}

	return max;

  }
 
 ///////////////////////////////////****** MAIN ******/////////////////////////////////////////


int main(){

//double alpha =pow(2.,-55);
  double alpha1=pow(2.,-55);
//4.26e-13;//pow(2.,-55);


  uint64_t n=1<<12;
  uint64_t pack=1; // data element size: 384 Byte
  uint64_t m=n/pack;
  
  uint64_t mm=1<<10; // CMUX per 2^10 data elts
  uint64_t tot =m/mm;// after 10 levels, run CMUX with tot(m/2^10) TRLWE samples

  cout << "the number of DB: "<< n << endl;
  cout << "the number of data in one plaintext(#packing): "<< pack << endl;
  cout << "the number of packed DB: "<< m  << endl;
  cout << "each size of data element: "<< "384Byte" << endl;
  cout << "small cmux gate with " << mm << "data elements" << endl;
  
 

 
  int index=rand()%m;
  cout << "-------------index:"<<index<<"------------------" << endl;
   //  rand()%m;//m-1;
  cout << "generating key switching materials (not depending on DB size) may take some time ... " << endl;


  Globals* env = new Globals();



	int32_t k = env->k;
	int32_t N = env->N;
	int32_t l=env->l;

  cout<< "the ciphertext polynomial degree: "<< N << endl;
    

  int32_t kpl =(k+1)*l;

  int Bgbit=env->bgbit;
  double log_2m=ceil(log(m)/log(2));
  Torus32 log2m = static_cast<int>(log_2m);
  double log_2mm=ceil(log(mm)/log(2)); 
  Torus32 log2mm = static_cast<int>(log_2mm);





  int MM1=12;
  uint64_t MM=int64_t(1)<<MM1; //plaintext modulus : 2^12

  double log_2t=ceil(log(MM)/log(2));
  uint64_t logt = static_cast<int>(log_2t);

  uint64_t db=logt*N*(m/64);
  
  uint64_t mdb = max(db,mm);
  uint64_t mct = max(mm,tot);



  TLweSample64 *enc_dataa =new TLweSample64(N);
  TLweSampleFFTa *enc_temp = new_array1<TLweSampleFFTa>(mct,N);
  TLweSampleFFTa *enc_data = new_array1<TLweSampleFFTa>(mct,N);
  










  Torus64Polynomial *data= new Torus64Polynomial(db);

  int64_t bit[log2m];

  int_to_bin_digit(index,log2m, bit);

	for(uint32_t j=0;j<db;j++)
	  data->coefs[j]=(int64_t(1)<<52)*(17*j*rand()%MM); // database as plaintext......?!?!?!?! with a little trick..........OMG!!!!


   Torus64Polynomial *indexdata= new Torus64Polynomial(N);
  
	  for(int32_t j=1;j<N;j++)
	  indexdata->coefs[j]=0;
	
	indexdata->coefs[0]=data->coefs[index];


//* client  is making a query.
//////////////////////////////////////////////////////////////////////////




  TGswSampleFFTa* extFFT= new_array1<TGswSampleFFTa>(log2m,l,N);
//  TGswSample64 *extract=new TGswSample64(l,N);
   TGswSample64 *extract = new_array1<TGswSample64>(log2m,l,N);

  LweSample64 *ext =new_array1<LweSample64>(log2m,N);
  LweSample64 *extbg =new_array1<LweSample64>(log2m,N);


  printf("Client gives a query\n\n");

  //auto start00 = std::chrono::high_resolution_clock::now();


  Torus64Polynomial *mess = new Torus64Polynomial(N);


  for (int i=0; i<N; i++) mess->coefs[i]=0;

    for (int32_t i = log2m-1; i>=0 ;--i)
    mess->coefs[log2m-1-i]=bit[i]*(UINT64_C(1) << (64-Bgbit));

  Torus64Polynomial *me = new Torus64Polynomial(N);
  for (int i=0; i<N; i++) me->coefs[i]=0;


    for (int32_t i = log2m-1; i>=0 ;--i)
        me->coefs[log2m-1-i]=bit[i]*(UINT64_C(1) << (64-2*Bgbit));

 TLweSample64 *bg =new TLweSample64(N);
 tLwe64Encrypt(bg, mess,alpha1, env);
 TLweSample64 *bgg =new TLweSample64(N);
 tLwe64Encrypt(bgg, me,alpha1, env);

	

	

	Torus64Polynomial *error_0 = new Torus64Polynomial(N);
	Torus64Polynomial *dec = new Torus64Polynomial(N);




 printf("log m ciphertexts(TGSWFFT form) are created.\n\n"); 
 printf("Server starts computation now\n\n");


  auto start0 = std::chrono::high_resolution_clock::now();
  for(int32_t i=0;i<log2m;i++){

	 //Sample Extraction

	tLweExtractLweSampleIndex64(&ext[i], bg, i, env);

	tLweExtractLweSampleIndex64(&extbg[i], bgg, i, env);
			}

  for(int s=0;s<log2m;++s){
    for (int j=0;j<2;++j){
	circuitPrivKS(&extract[s].allsamples[2*j],j, &ext[s],env);
        circuitPrivKS(&extract[s].allsamples[2*j+1],j, &extbg[s], env);

         }
	}
  auto end0 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time0 = end0-start0;
  std::cout << "Query unpacking step takes : "<<execution_time0.count()<< " ms " <<std::endl;

  
     
////////////////////////////////////////////////////////////////////////////////



  printf("STep2 starts\n\n");
  



  auto start = std::chrono::high_resolution_clock::now();

  auto start20 = std::chrono::high_resolution_clock::now();
  for (int s=0;s<log2m;s++){
     	 for (int i=0;i<kpl;i++)
   		for (int q=0;q<=k;q++)
 		     TorusPolynomial64_ifft_lvl2(&extFFT[s].allsamples[i].a[q],&extract[s].allsamples[i].a[q],  env);
            }

  auto end20 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time20 = end20-start20;
  std::cout <<"query:conversion to FFT: " << execution_time20.count()<<" ms "<<std::endl;

  int32_t j=0;
   int32_t z =0;
   int32_t mk=mm;
 
  while(z<tot){

	  while(j<mm)
 	 {
 		  
 	 CMuxFFTdb(&enc_temp[j],&extFFT[0],data->coefs[z*mk+j],data->coefs[z*mk+j+1], env);
// CMuxFFT(&enc_temp[j],&extFFT[0],&enc_data[j],&enc_data[j+1], env);

	
	j=j+2;
       } 	  
 	//printf("%d\n",z*mk);	 

	j=0;

	mm=mm/2;
	

	for (int32_t i=1; i<log2mm;++i)
 	{

       

 	  while(j<mm)
  	   {
  	   CMuxFFTa(&enc_temp[j],&extFFT[i],&enc_temp[j*2],&enc_temp[(j+1)*2], env);
	   

  	    j=j+2;



  	   }
	

  	j=0;
 	mm=mm/2;

  	}


 
   for(int q=0;q<=k;++q)
  
    {for(int c=0;c<N;++c)	
     enc_data[z].a[q].values[c]=enc_temp[0].a[q].values[c];
    }
    mm=mk;
    j=0;
    z=z+1;
	



  }


  
j=0;



 while(j<tot)
 	{
 		  
 //	 CMuxDecompFFTa(&enc_temp[j],&extFFT[0],decompFFT+(kpl*j), env);
	 CMuxFFTa(&enc_temp[j],&extFFT[log2mm],&enc_data[j],&enc_data[j+1], env);


 	   j=j+2;

	  }
 		 
	j=0;
	tot=tot/2;
	

	for (int32_t i=log2mm+1; i<log2m;++i)
 	{

       

 	  while(j<tot)
  	   {
  	   CMuxFFTa(&enc_temp[j],&extFFT[i],&enc_temp[j*2],&enc_temp[(j+1)*2], env);


  	    j=j+2;



  	   }
	

  	j=0;
 	tot=tot/2;

 }






  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time = end-start;
  std::cout <<"main computation step(running m-1 cmux gates) takes: " << execution_time.count()<<" ms "<<std::endl;


   

//After the second step the output ciphertext is packed including the desired item.








  for (int q = 0; q <= k; ++q) 	

		TorusPolynomial64_fft_lvl2(&enc_dataa->a[q],&enc_temp[0].a[q], env);

   Torus64Polynomial *decrypt = new Torus64Polynomial(N);
  auto start11 = std::chrono::high_resolution_clock::now();
 

	tLwe64Phase_lvl2(decrypt,enc_dataa, env);


  auto end11 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double,std::milli> execution_time11 = end11-start11;
  std::cout <<" answer decode takes: "<< execution_time11.count()<<" ms "<<std::endl;
	


  Torus64Polynomial *error = new Torus64Polynomial(N);
   
   for (int i=0;i<N;++i){

	error->coefs[i]=decrypt->coefs[i]-indexdata->coefs[i];
    }
  int64_t temp1 =  array_max(error,N);
 	//infinite norm of add error
   	
  double bit_ea = ceil(log(temp1)/log(2));
  cout<< "output noise budget: " << 64-logt-bit_ea-1 << endl;




  int aa=0;
  cout << "Test TLweSymDecrypt " << endl;
    for (int32_t i = 0; i < N; ++i) {
         if (abs(decrypt->coefs[i]-indexdata->coefs[i]) >int64_t(1)<<(63-MM1))

   aa+=1;
    }
	if (aa>0)
    cout<< "decryption failure? = "<< "Yes" << endl;
	else
    cout<< "decryption failure? = "<< "No" << endl;
    






}
