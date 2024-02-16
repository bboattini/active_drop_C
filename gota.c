/*********************************************************************************************
*																							 *
*        Programa Potts Heitor modificado para incluir a mistura água/óleo                   *
*			        Modificado por Cristina Gavazoni, Junho 2020							 *
*																							 *
*********************************************************************************************/
/* Não está implementado a continuação de runs antigos!!!*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mc.h"
#include "pointers.h"
#ifdef USEGFX
#include "board3d.h"
#endif

/*********************************************************************************************
*                                 Declaração das funções                                     *
*********************************************************************************************/

void create_pointers(void);
void openfiles(void);
void initialization(void);
void initial_state(int *s, int l, int rg, int initialstate);
void xyz(int *s, int l, int m);
void dynamics (int *s,int num_steps, double Aw);
void flip_spin (long site, long new_s);
void add_to_interface(int site);
void remove_from_interface(int site);
void free_pointers(void);
void measure_angle(int num_steps);
double calculate_energy(void);
void save_conf(int num_steps, int iout);
void teste_flip_droplet(int site,int s_teste,double delta_s, double delta);
void teste_flip_surface(int site,int s_teste,double delta_s, double delta);


/*********************************************************************************************
*                       Declarando parâmetros da simulação - técnicos                        *
*********************************************************************************************/

#define mc_steps   	       1000  // Número de passos de MC totais
#define n_mesure       	   10   // Intervalo para salvar medidas
#define n_teste       	   99990   // Intervalo para salvar medidas

#define temp           	   13.0  // Temperatura
#define kB             	   1.0  // Constante de Boltzman

#define t_neigh        	   26   // Número de vizinhos
#define t_close_neigh  	   18   // Número de vizinhos próximos 

#define NUM_CONF            1 // DEIXAR IGUAL A 1 !!

/*********************************************************************************************
*                    Declarando Parâmetros constantes durante a simulação                    *
*********************************************************************************************/

#define frac_h 				1.0 // Se > 1 a gota flutua sobre a superfície
#define h_base       		3 //Altura da Base

#define GW                  0.0001 // Massa agua em um sítio
#define GO                  0.000077 // Mass óleo em um sítio
#define Lambda_w            0.01 // Parâmetro de volume
#define Lambda_o            0.01 // Parâmetro de volume

#define eps_SG              0.96 // E superficial S-G

#define eps_WG              2.70 // E superficial L-G agua
#define cos_theta_w         (cos(111.0*M_PI/180.0))
#define eps_SW              (eps_SG-eps_WG*cos_theta_w)// E superficial S-L agua


#define eps_OG              1.04 // E superficial L-G oleo
#define cos_theta_O         (cos(53.0*M_PI/180.0))
#define eps_SO              (eps_SG-eps_OG*cos_theta_O)// E superficial S-L oleo

#define eps_WO              2.06 // E superficial óleo-água

/*********************************************************************************************
*                                   Declarando das variáveis                                 *
*********************************************************************************************/

unsigned long identifier;

//--------------------------------------------------------------------------------------------
// Variáveis da configuração inicial da gota
//--------------------------------------------------------------------------------------------

int l, l2, l3, lm;  // Tamanho do sistema e multiplos
double r0;
int rg, rg2, rg3;   // Tamanho da gota e multiplos
int h, w,a, d;     	// Altura, largura, separaco pilar e d=w+a
int initialstate;  	// Estado inicial: W ou CB
double hm, thmed;  	// h médio dos pilares, ângulo médio
double fvmed_w,fvmed_o; // fração med gota de agua e óleo
int c_theta, c_fv_w, c_fv_o;// contador angulo, contador das porcentagens de volume
int hmax, hmin;    	// multiplicador para numero de blocos, h máximo e mínimo da conf

double v0, v0_w, v0_o; // volumes desejados
int t_vol, t_vol_w, t_vol_o; // volumes desejados
double f_w, f_o;                // Frações de água e óleo
double fw, fo;                // Frações de água e óleo

double T;

char CI[4]; // Guarda o nome do arquivo

//--------------------------------------------------------------------------------------------
// Ponteiros
//--------------------------------------------------------------------------------------------

int *s ; // Spins
int *nv_cw,*nv_co, *nv_w, *nv_o, *nv_s, *nv_g; // Números de vizinhos de cada tipo
int *inter_pos, *w_inter, *inter_mtx; // Sítios na interface

//--------------------------------------------------------------------------------------------
// Variaveis usadas na simulação
//--------------------------------------------------------------------------------------------

int n_w, n_o, n_s; //número de sítios água e óleo
int int_label; //label da lista de sítios na interface
int vol, vol_w, vol_o; // volumes calculados
double Gw, Go, Aw;

//--------------------------------------------------------------------------------------------
// Variaveis dos observáveis
//--------------------------------------------------------------------------------------------

double angle_x,angle_ellipse_x,base_radius_x,drop_radius_x,major_axis_x,minor_axis_x;
double angle_y,angle_ellipse_y,base_radius_y,drop_radius_y,major_axis_y,minor_axis_y;
int centre_x,centre_y,drop_height,drop_bottom;
int num_grooves,num_wenzel;

double angle_x_o,angle_ellipse_x_o,base_radius_x_o,drop_radius_x_o,major_axis_x_o, minor_axis_x_o;
double angle_y_o,angle_ellipse_y_o,base_radius_y_o,drop_radius_y_o,major_axis_y_o, minor_axis_y_o;
int centre_x_o,centre_y_o,drop_height_o,drop_bottom_o;
int num_grooves_o,num_wenzel_o;

//--------------------------------------------------------------------------------------------
// Nomes dos arquivos de entrada e saida 
//--------------------------------------------------------------------------------------------

char output_file1[100],input_file[100]; //Nome dos arquivos
FILE *fconf,*fp1,*finit,*fxyz,*flast,*fp_in,*fout,*foil,*fwat; // Localizadores dos arquivos

//--------------------------------------------------------------------------------------------
// Estatística
//--------------------------------------------------------------------------------------------

long int trial,accepted,acc_interface,acc_teste;

//--------------------------------------------------------------------------------------------
// Coisas para multiplas simulações
//--------------------------------------------------------------------------------------------

int cont; //controla o numero de continuacoes de uma dada simualcao
int file_in=0,confstep,iconf;
long int tempo_in=0;
unsigned int numsteps;

/*********************************************************************************************
==============================================================================================
=                                      PROGRAMA PRINCIPAL                                    =
==============================================================================================
*********************************************************************************************/
int main(int argc,char *argv[])
{

	int i;

//--------------------------------------------------------------------------------------------
//  Lendo os parâmetros de entrada e definindo parametros da conf inicial
//--------------------------------------------------------------------------------------------

	identifier = 0;

	if(argc==17) 
	{
   		 for (i=1;i<argc;i++) 
		{
      		if (!strcmp(argv[i],"-L"))       l=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-R"))  r0=atof(argv[++i]);
      		else if (!strcmp(argv[i],"-a"))  a=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-h"))  h=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-w"))  w=atoi(argv[++i]);
			// Substituir Aw por f_o para recuperar a forma original
      		else if (!strcmp(argv[i],"-fo")) Aw=atof(argv[++i]);
      		else if (!strcmp(argv[i],"-CI")) initialstate=atoi(argv[++i]);
      		else if (!strcmp(argv[i],"-s"))  seed=atol(argv[++i]);
      		else {
				fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",argv[i]);
				exit(-1);
      			}/* Encerra Else */
    	} /* Encerra FOR */

    	if(initialstate==1) {
    	  sprintf(CI,"CB_");
    	}
    	else if(initialstate==2) {
    	  sprintf(CI,"WE_");
    	}
    	else {
		exit(-1);
    	}

		rg = (int)r0;
		f_o = 0.0;

		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f.out",CI,l,rg,a,h,w,Aw);	
		fout = fopen(output_file1,"w");	

		fflush(fout);

/*		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_oil.out",CI,l,rg,a,h,w,Aw);	*/
/*		foil = fopen(output_file1,"w");	*/
/*		fprintf(foil ,"#  site   s[site]    s_teste    nw     no    ng    ns      ei     ef     ef-ei   delta_s     delta\n");*/

/*		fflush(foil);*/

/*		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_water.out",CI,l,rg,a,h,w,Aw);	*/
/*		fwat = fopen(output_file1,"w");	*/
/*		fprintf(fwat ,"#  site   s[site]    s_teste    nw     no    ng    ns      ei     ef     ef-ei   delta_s     delta\n");*/

/*		fflush(fwat);*/


		rg2 = rg*rg;
		rg3 = rg2*rg;
	
		l2 = l*l;
		l3 = l2*l;
		lm = l/2; 

		f_w =  1.0 - f_o;

		v0 = (4.0*M_PI*r0*r0*r0/3.0);
		v0_w = f_w*v0;
		v0_o = f_o*v0;

		t_vol = (int) v0;
		t_vol_w = (int) v0_w;
		t_vol_o = (int) v0_o;

		T = temp;

		fprintf(fout,"=====================================================================\n");
		fprintf(fout,"=                 Parâmetros Iniciais do Sistema                    =\n");
		fprintf(fout,"=====================================================================\n");
		fprintf(fout,"L      : %d\n",l);
		fprintf(fout,"R      : %d\n",rg);
		fprintf(fout,"V      : %d\n",t_vol);
    	fprintf(fout,"V_w    : %d\n",t_vol_w);
    	fprintf(fout,"V_o    : %d\n",t_vol_o);
    	fprintf(fout,"fw     : %f\n",f_w);
		// Substituir Aw por f_o para recuperar a forma original
    	fprintf(fout,"fo     : %f\n",Aw);
    	fprintf(fout,"h      : %d\n",h);
    	fprintf(fout,"w      : %d\n",w);
    	fprintf(fout,"a      : %d\n",a);
    	fprintf(fout,"State  : %d\n",initialstate);
    	fprintf(fout,"eps_SG : %f\n", eps_SG);
    	fprintf(fout,"eps_WO : %f\n", eps_WO);
    	fprintf(fout,"eps_SW : %f\n", eps_SW);
    	fprintf(fout,"eps_WG : %f\n", eps_WG);
    	fprintf(fout,"eps_SO : %f\n", eps_SO);
    	fprintf(fout,"eps_OG : %f\n", eps_OG);	


    	fprintf(fout,"CI     : %s\n",CI);
    	fprintf(fout,"Seed   : %ld\n",seed);
		fprintf(fout,"=====================================================================\n");
		fprintf(fout,"\n");

	} /* Encerra IF */

	else {
    	fprintf(stderr,"Error.     Invalid outout. Use:\n");
    	fprintf(stderr,"             -a interpillar distance\n");
    	fprintf(stderr,"             -h pillar height\n");
    	fprintf(stderr,"             -L l\n");
    	fprintf(stderr,"             -R LWATER\n");
    	fprintf(stderr,"(optional)   -f ID from file to continue simulation\n");
    	fprintf(stderr,"(optional)      BE CAREFULL WITH INPUT PARAMETERS!\n");
    	exit(-1);
  	} //Fim o ELSE

//--------------------------------------------------------------------------------------------
//  Escolhendo uma semente aleatóriamente
//--------------------------------------------------------------------------------------------

	identifier = start_randomic_seed_ext(seed);
	if(identifier!=seed) {printf("id: %ld s: %ld\n",identifier,seed);exit(-1);}

//--------------------------------------------------------------------------------------------
//  Criando Ponteiros
//--------------------------------------------------------------------------------------------

	create_pointers();

//--------------------------------------------------------------------------------------------
//  Abrindo arquivos
//--------------------------------------------------------------------------------------------

	openfiles();

//--------------------------------------------------------------------------------------------
//  Gerando a Configuração inicial
//--------------------------------------------------------------------------------------------

	initialization();
	xyz(s,l,0);
	save_conf(0,0);

//--------------------------------------------------------------------------------------------
//  Estatística
//--------------------------------------------------------------------------------------------

	trial=0; accepted=0;acc_interface=0;acc_teste=0;

//--------------------------------------------------------------------------------------------
//  Dinâmica e medida dos observáveis
//--------------------------------------------------------------------------------------------


	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                       Energias do Sistema                         =\n");
	fprintf(fout,"=====================================================================\n");

	thmed = 0;
	c_theta=0;

	fvmed_w = 0;
	fvmed_o = 0;
	c_fv_w = 0;
	c_fv_o = 0;	

	for (i=0; i<=mc_steps; ++i) 
	{
		Gw=GW;
		Go=GO;


		dynamics(s,i,Aw);

		fw = (float)vol_w/(vol);
		fo = (float)vol_o/(vol);

		if(i%n_mesure==0) 
		{
			measure_angle(i);
			save_conf(i,0);
			
		}//fim do IF

	}//Fim do FOR

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n");

	save_conf(mc_steps,1);
	xyz(s,l,1);

//--------------------------------------------------------------------------------------------
//  Estatística
//--------------------------------------------------------------------------------------------

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                         Fim da simulação                          =\n");
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"Número de partículas de água    : %d\n",vol_w);
	fprintf(fout,"Número de partículas de óleo    : %d\n",vol_o);
	fprintf(fout,"Número de partículas de líquido : %d\n",vol);
	fprintf(fout,"Número de sítios na interface   : %d\n",int_label);
	fprintf(fout,"Fração de óleo                  : %f\n",fo);
	fprintf(fout,"Fração de água                  : %f\n",fw);
	fprintf(fout,"Ângulo de contato			    : %f\n",thmed/c_theta);
	fprintf(fout,"Fração de água na gota		    : %f\n",fvmed_w/c_fv_w);
	fprintf(fout,"Fração de óleo na gota		    : %f\n",fvmed_o/c_fv_o);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n");
	fflush(stdout);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                            Estatísticas                           =\n");
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"Trial: %ld\n",trial);
	fprintf(fout,"Accepted: %ld     Ratio: %f\n",accepted,(float)accepted/trial);
	fprintf(fout,"Interface: %ld    Ratio: %f\n",acc_interface,(float)acc_interface/trial);
	fprintf(fout,"Teste: %ld        Ratio: %f\n",acc_teste,(float)acc_teste/trial);
	fprintf(fout,"Variação do volume: %f\n", (float)vol/t_vol);
	fprintf(fout,"Variação fo: %f   Variação fw: %f\n",fo/f_o, fw/f_w);
	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"\n\n\n");

//--------------------------------------------------------------------------------------------
//  Encerrando
//--------------------------------------------------------------------------------------------

 	free_pointers();
 	fclose(fp1);

} //Encerra o Programa
/********************************************************************************************/



/*********************************************************************************************
==============================================================================================
=                                    Subrotinas e Funções                                    =
==============================================================================================
*********************************************************************************************/

/*********************************************************************************************
*                             Rotina que abre os arquivos de saída                            *
*********************************************************************************************/
void openfiles(void)
{
	char lixo;
	
	file_in = 1;

	// Nome do arquivo - o identifier tem que estar correto!
	sprintf(input_file,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_LAST.dsf",CI,l,rg,a,h,w,Aw);
  
 	 // Abrindo entrada
	if((fp_in = fopen(input_file,"r")) == NULL)
	{
      file_in = 0;
    } // Fim do IF

//--------------------------------------------------------------------------------------------
//   Arquivos neceessários para continuar um run anterior - Olhar melhor
//--------------------------------------------------------------------------------------------
	else 
	{
		if((fscanf(fp_in,"%c %ld %d",&lixo,&tempo_in,&confstep)!=3)) 
		{
		fprintf(stderr,"\n\n ERRO lendo arq de entrada: primeira linha\n\n");
		exit(-1); 
    	}// Fim do IF


    	fprintf(stderr,"Último tempo: %ld\n",tempo_in);
    	fprintf(stderr,"Última amostra: %d\n",confstep);
    	fprintf(stderr,"Lixo: %c\n",lixo);
    	numsteps = tempo_in;
    	if(numsteps==mc_steps) 
		{
			fprintf(stderr,"\n\nESTA AMOSTRA JA ACABOU! SE QUISER AUMENTAR O TEMPO DE \n");
			fprintf(stderr,"SIMULACAO, AUMENTE A VARIAVEL TOTALSTEPS!!! \n\n");
			exit(1);
		}// Fim do IF

	    numsteps++;
	    fclose(fp_in);
	    // leio os sitios da interface
	    fflush(stdout); 
  
    	// Nome do arquivo de saída 
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f.dsf",CI,l,rg,a,h,w,Aw);
    
		fp1 = fopen(output_file1,"a");
		fprintf(fp1,"### Continuando..\n");
	    fflush(fp1);
    
	    fflush(stdout);
   
	}// Fim do ELSE

//--------------------------------------------------------------------------------------------
//   Arquivos novos - Para NOVOS runs
//--------------------------------------------------------------------------------------------

	if(file_in==0)  // Se é uma simulação nova
	{
		// Arquivo dsf - Observáveis  
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f.dsf",CI,l,rg,a,h,w,Aw);

		fflush(stdout);
		fp1 = fopen(output_file1,"w"); //fp1 é o localizador desse arquivo

		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"# =                     Parâmetros da Simulação                       =\n");
		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"# Gota - CB (1) ou Wenzel(2)  3D : %d\n",initialstate);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# L = %d   Rg = %d   fo = %3.2f\n" ,l,rg, f_o);
		fprintf(fp1,"# w = %d      h = %d    a = %d\n",w,h,a);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Gw = %5.4f   Lambda_w = %5.4d\n",Gw, Lambda_w);
		fprintf(fp1,"# Go = %5.4f   Lambda_o = %5.4d\n",Go, Lambda_o);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# eps_SG = %3.2f     eps_WO = %3.2f\n",eps_SG, eps_WO);
		fprintf(fp1,"# eps_SW = %3.2f     eps_WG = %3.2f\n",eps_SW, eps_WG);
		fprintf(fp1,"# eps_SO = %3.2f     eps_OG = %3.2f\n",eps_SO, eps_OG);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Temp = %4.2f           kB = %3.2f\n",temp,kB);
		fprintf(fp1,"# \n");
		fprintf(fp1,"# Total time = %d\n",mc_steps);
		fprintf(fp1,"# =====================================================================\n");
		fprintf(fp1,"#  t        V      Vw    Vo       E           theta_x  theta_y     R_x        R_y    B px  py   B_x         B_y     H   vb\n");
   		fflush(fp1);
    
		fflush(stdout);

      // Arquivo de Conficguração
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_conf.dsf", CI,l,rg,a,h,w,Aw);

		fconf = fopen(output_file1,"w");

		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# =                     Parâmetros da Simulação                       =\n");
		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# Gota - CB (1) ou Wenzel(2)  3D : %d\n",initialstate);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# L = %d   Rg = %d   fo = %3.2f\n" ,l,rg, f_o);
		fprintf(fconf,"# w = %d      h = %d    a = %d\n",w,h,a);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Gw = %5.4f   Lambda_w = %5.4d\n",Gw, Lambda_w);
		fprintf(fconf,"# Go = %5.4f   Lambda_o = %5.4d\n",Go, Lambda_o);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# eps_SG = %3.2f     eps_WO = %3.2f\n",eps_SG, eps_WO);
		fprintf(fconf,"# eps_SW = %3.2f     eps_WG = %3.2f\n",eps_SW, eps_WG);
		fprintf(fconf,"# eps_SO = %3.2f     eps_OG = %3.2f\n",eps_SO, eps_OG);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Temp = %4.2f           kB = %3.2f\n",temp,kB);
		fprintf(fconf,"#\n");
		fprintf(fconf,"# Total time = %d\n",mc_steps);
		fprintf(fconf,"# =====================================================================\n");
		fprintf(fconf,"# t and interface's sites\n");

		fflush(fconf);

		fflush(stdout);

/*#endif*/

		// Arquivo xyz init
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_init.xyz",CI,l,rg,a,h,w,Aw);	
		finit = fopen(output_file1,"w");	

		fflush(finit);

		// Arquivo xyz
        sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f.xyz",CI,l,rg,a,h,w,Aw);	
		fxyz = fopen(output_file1,"w");	

		fflush(fxyz);
        

  } // Fim do IF

  return;
}
/*********************************************************************************************
*                              Rotina que inicializa o programa                              *
*********************************************************************************************/
void initialization(void)
{
	int i;
  
	initial_state(s,l,rg,initialstate);
  
	n_w=0;
	n_o=0;
	n_s=0;
	int_label=0;
	
	for (i=0; i < l3; ++i) 
    {
		if (s[i]==1) ++n_w;
		if (s[i]==2) ++n_o;
		if (s[i]==9) ++n_s;

		if (inter_pos[i]==1) 
		{
			w_inter[int_label]=i;
			inter_mtx[i]=int_label;
			++int_label;

		} // Fim do If

		else inter_mtx[i]=-1;

    } //Fim do FOR

	vol_w = n_w;
	vol_o = n_o;
	vol = vol_w + vol_o;

	fprintf(fout,"=====================================================================\n");
	fprintf(fout,"=                        Início da simulação                        =\n");
	fprintf(fout,"=====================================================================\n");   
	fprintf(fout,"Número de partículas de água    : %d\n",vol_w);
	fprintf(fout,"Número de partículas de óleo    : %d\n",vol_o);
	fprintf(fout,"Número de partículas de líquido : %d\n",vol);
	fprintf(fout,"Número de sítios na interface   : %d\n",int_label);
	fprintf(fout,"H médio                         : %f\n",hm);
	fprintf(fout,"H máximo                        : %d\n",hmax - h_base);
	fprintf(fout,"H mínimo                        : %d\n",hmin - h_base);
	fprintf(fout,"=====================================================================\n"); 
	fprintf(fout,"\n"); 
	fflush(stdout);

  return;
}
/*********************************************************************************************
*                          Rotina gera a gota com a distribuição desejada                    *
*********************************************************************************************/
void initial_state(int *s, int L, int Rg, int Initialstate)  
{
  
	int i,j,k;
	int l2,l3,rg2;
	int site;
 	long int count=0;
	int neigh[26];
	int x,y,z,xn,yn,zn;
	int zm, z_teste, nm;
	double fator;
/*	int nw, nt,*/
	double teste;

	l2=L*L;
	l3=l2*L;

	rg2=Rg*Rg; 

/*	nw = t_vol_w -1; // usado pra preencher a gota modo ordenado*/

	for (i=0;i<l3;++i)
    {
      s[i]=0;
      nv_cw[i]=0; //vizinhos proximos água
      nv_co[i]=0; //vizinhos proximos oleo
      nv_s[i]=0; //vizinhos solido
      nv_w[i]=0; //vizinhos agua
      nv_o[i]=0; //vizinhos oleo
      nv_g[i]=0; //vizinhos gas
      inter_pos[i]=0; 
      w_inter[i]=-1;

    } //Fim do FOR

//--------------------------------------------------------------------------------------------
//   Criando uma nova gota
//--------------------------------------------------------------------------------------------
	if(file_in==0) 
	{

		numsteps=0;
		confstep=0;
		count = 0;
/*		nt=0;*/

		switch (initialstate)
		{
      	case 1 : // CB

			for (k = 0; k < L; ++k)
				for(j = 0; j < L; ++j)
					for(i = 0; i < L; ++i)
			{

					long z=frac_h*Rg+h+h_base+1;
					site=k*l2+L*j+i;
					if ( (i-L/2+1)*(i-L/2+1) + (j-L/2+1)*(j-L/2+1)+(k-z+1)*(k-z+1) <= rg2)
		  			{

						if ( f_o == 0.0) //Se for apenas água
		  				{
		    				*(s+site)=1;
		    				count++;
						}//fim do IF
						else if ( f_o == 1.0)//Se for apenas óleo
		  				{
		    				*(s+site)=2;
		    				count++;
						}//fim do ELSE IF

						else //se for uma mistura
						{
							/*---------------------------------*/ 
							/*Para gerar uma conf aleatória*/ 

							teste = FRANDOM;
							if ( (teste>f_o)  )
							{

								*(s+site)=1;
		    					count++;

							} //fim do IF
							else if ( (teste<=f_o)  )
							{

		    					*(s+site)=2;
		    					count++;

							} //fim do ELSE IF
							/*---------------------------------*/

							/*-------------------------------*/
							/*Para gerar uma conf ordenada   */
/*							if (nt < nw)*/
/*							{*/
/*		    					*(s+site)=1;*/
/*								nt++;*/
/*		    					count++;*/
/*							} //fim do IF	*/
/*							else*/
/*							{*/
/*		    					*(s+site)=2;*/
/*								nt++;*/
/*		    					count++;						*/
/*							}//fim do ELSE*/
							/*-------------------------------*/
							
						}//fim do ELSE

		  			} //fim do IF

			} //Fim do FOR
		break;

		case 2 : // We

			fator = pow(2.0,0.666);
			for (k = 0; k < L; ++k)
				for(j = 0; j < L; ++j)
					for(i = 0; i < L; ++i)
			{

					long z=h_base;
					site=k*l2+l*j+i;
					if ( (i-L/2+1)*(i-L/2+1) + (j-L/2+1)*(j-L/2+1)+(k-z+1)*(k-z+1) <= fator*rg2)
		  			{
		    			if ( f_o == 0.0) //Se for apenas água
		  				{
		    				*(s+site)=1;
		    				count++;
						}//fim do IF
						else if ( f_o == 1.0)//Se for apenas óleo
		  				{
		    				*(s+site)=2;
		    				count++;
						}//fim do ELSE IF

						else //se for uma mistura
						{
							/*---------------------------------*/ 
							/*Para gerar uma conf aleatória*/ 

							teste = FRANDOM;
							if ( (teste>f_o)  )
							{

								*(s+site)=1;
		    					count++;

							} //fim do IF
							else if ( (teste<=f_o)  )
							{

		    					*(s+site)=2;
		    					count++;

							} //fim do ELSE IF
							/*---------------------------------*/

							/*-------------------------------*/
							/*Para gerar uma conf ordenada   */
/*							if (nt < nw)*/
/*							{*/
/*		    					*(s+site)=1;*/
/*								nt++;*/
/*		    					count++;*/
/*							} //fim do IF	*/
/*							else*/
/*							{*/
/*		    					*(s+site)=2;*/
/*								nt++;*/
/*		    					count++;						*/
/*							}//fim do ELSE*/
							/*-------------------------------*/
							
						}//fim do ELSE
		  			} //fim do IF

			}//fim do FOR
		break;

		}//Fim do SWITCH

	}//Fim do IF

//--------------------------------------------------------------------------------------------
//   Criando a superfície de pilares
//--------------------------------------------------------------------------------------------

	for (k=0;k<h+h_base;++k) 
	{
		for (j=0;j<L;++j) 
		{
			for (i=0;i<L;++i) 
			{ 
				
				site=k*l2+j*L+i;

				// Superfície de pilares
				if ( (k<h_base) || ((i%(w+a)<w)&&(j%(w+a)<w)) ) s[site]=9;
	
      		} //fim do FOR
		}//fim do FOR   
	} //fim do FOR

//--------------------------------------------------------------------------------------------
// Calculando h médio e hmax
//--------------------------------------------------------------------------------------------

	zm = 0;
	nm = 0;
	hmax = 0;
	hmin = h_base+1;

	for (i=0;i<L;++i) 
		{
			for (j=0;j<L;++j) 
			{

				z_teste = 0;
				for (k=0;k<L;++k) 
				{ 

					site=k*l2+L*j+i;
					if((s[site]==9) && (k>z_teste)) z_teste=k;
					if((s[site]==9) && (k>hmax)) hmax=k;
	
      			} //fim do FOR
				
				zm = zm + z_teste;
				nm++;

			}//fim do FOR   
		} //fim do FOR  

	zm = zm - h_base;
	hm = (float)zm/nm;  

//--------------------------------------------------------------------------------------------
//   Criando a lista de vizinhos
//--------------------------------------------------------------------------------------------	   
	for(i=0;i<l3;++i)
    {
		//vizinhos
		site = i;
		x = site%L;
		y = (site/L)%L;
		z = site/l2;
      
		// 0
		yn = (y==0)?L-1:y-1;
		neigh[0]=x + yn*L + z*l2;
		// 1
		yn = (y==L-1)?0:y+1;
		neigh[1] = x + yn*L + z*l2;
		// 2
		xn = (x==L-1)?0:x+1;
		neigh[2] = xn + y*L + z*l2;
		// 3
		xn = (x==0)?L-1:x-1;
		neigh[3] = xn + y*L + z*l2;
		// 4
		zn = (z==0)?L-1:z-1;
		neigh[4] = x + y*L + zn*l2;
		// 5
		zn = (z==L-1)?0:z+1;
		neigh[5] = x + y*L + zn*l2;
		//6 
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		neigh[6] = xn + yn*L + z*l2;
		// 7
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		neigh[7] = xn + yn*L + z*l2;
		// 8 
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		neigh[8] = xn + yn*L + z*l2;
		// 9
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		neigh[9] = xn + yn*L + z*l2;
		// 10
		xn = (x==L-1)?0:x+1;
		zn = (z==0)?L-1:z-1;
		neigh[10] = xn + y*L + zn*l2;
		// 11
		xn = (x==0)?L-1:x-1;
		zn = (z==0)?L-1:z-1;
		neigh[11] = xn + y*L + zn*l2;
		// 12
		xn = (x==L-1)?0:x+1;
		zn = (z==L-1)?0:z+1;
		neigh[12] = xn + y*L + zn*l2;
		// 13
		xn = (x==0)?L-1:x-1;
		zn = (z==L-1)?0:z+1;
		neigh[13] = xn + y*L + zn*l2;
		// 14
		zn = (z==0)?L-1:z-1;
		yn = (y==0)?L-1:y-1;
		neigh[14] = x + yn*L + zn*l2;
		// 15
		zn = (z==L-1)?0:z+1;
		yn = (y==0)?L-1:y-1;
		neigh[15] = x + yn*L + zn*l2;
		// 16
		zn = (z==0)?L-1:z-1;
		yn = (y==L-1)?0:y+1;
		neigh[16] = x + yn*L + zn*l2;
		// 17
		zn = (z==L-1)?0:z+1;
		yn = (y==L-1)?0:y+1;
		neigh[17] = x + yn*L + zn*l2;
		// 18
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		zn = (z==0)?L-1:z-1;
		neigh[18] = xn + yn*L + zn*l2;
		// 19
		xn = (x==L-1)?0:x+1;
		yn = (y==0)?L-1:y-1;
		zn = (z==L-1)?0:z+1;
		neigh[19] = xn + yn*L + zn*l2;
		// 20
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		zn = (z==0)?L-1:z-1;
		neigh[20] = xn + yn*L + zn*l2;
		// 21
		xn = (x==0)?L-1:x-1;
		yn = (y==0)?L-1:y-1;
		zn = (z==L-1)?0:z+1;
		neigh[21] = xn + yn*L + zn*l2;
		// 22
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		zn = (z==0)?L-1:z-1;
		neigh[22] = xn + yn*L + zn*l2;
		// 23
		xn = (x==L-1)?0:x+1;
		yn = (y==L-1)?0:y+1;
		zn = (z==L-1)?0:z+1;
		neigh[23] = xn + yn*L + zn*l2;
		// 24
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		zn = (z==0)?L-1:z-1;
		neigh[24] = xn + yn*L + zn*l2;
		// 25	  
		xn = (x==0)?L-1:x-1;
		yn = (y==L-1)?0:y+1;
		zn = (z==L-1)?0:z+1;
		neigh[25] = xn + yn*L + zn*l2;
 
//--------------------------------------------------------------------------------------------
//   Contando o numero de vizinhos de cada tipo
//--------------------------------------------------------------------------------------------	
		for (j=0;j<t_close_neigh;++j)
		{

			if (s[neigh[j]]==1) 
	    	{
				++nv_cw[i];
				++nv_w[i];
	    	} //fim do IF

			else if (s[neigh[j]]==2) 
	    	{
				++nv_co[i];
				++nv_o[i];
	    	} //fim do ELSE IF

			else if (s[neigh[j]]==0) 
	    	{
				++nv_g[i];
	    	}//fim do ELSE IF

			else if (s[neigh[j]]==9)
	    	{
				++nv_s[i];
	    	}//fim do ELSE IF

		}//fim do FOR

//-----------------------------------------------------------------
		for (j=t_close_neigh;j<t_neigh ;++j)
		{
			if (s[neigh[j]]==1) 
	    	{
				++nv_w[i];
	    	} //fim do IF

			else if (s[neigh[j]]==2) 
	    	{
				++nv_o[i];
	    	} //fim do ELSE IF

			else if (s[neigh[j]]==0) 
	    	{
				++nv_g[i];
	    	}//fim do ELSE IF

			else if (s[neigh[j]]==9)
	    	{
				++nv_s[i];
	    	}//fim do ELSE 

		}//fim do FOR

//--------------------------------------------------------------------------------------------
//   Vendo se o sitio esta na interface
//--------------------------------------------------------------------------------------------	
		for (j=0;j<t_close_neigh;++j)
		{

			if ( (s[i]+s[neigh[j]] == 1)  ) //Gas + Água
		    {
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do IF
			else if ( (s[i]+s[neigh[j]] == 2) && (s[i]!=s[neigh[j]]) ) //Gas + óleo
		    {	
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do ELSE IF
			else if ( (s[i]+s[neigh[j]] == 3) ) //Água + óleo
		    {
				inter_pos[i]=1;
				inter_pos[neigh[j]]=1;
			}//fim do ELSE IF
			else if ( (s[i] == 1)  && (s[neigh[j]] == 9) )  //Água + solido
			{  
				inter_pos[i]=1;
	  		}//fim do ELSE IF
			else if ( (s[i] == 2)  && (s[neigh[j]] == 9)  )  //Oleo + solido
			{  
				inter_pos[i]=1;
	  		}//fim do ELSE IF
	 
		}//fim do FOR
//--------------------------------------------------------------------------------------------

    }//fim do FOR

 return;

}
/*********************************************************************************************
*                              Rotina que escreve o arquivo xyz                              *
*********************************************************************************************/
void xyz(int *s, int l, int m)  
{

	int i;
	int l2,l3; 
	int site, npart;
	int x,y,z;

	l2 = l*l;
	l3 = l2*l;

	npart = n_w+n_o+n_s;

	if ( (m == 0) ) 
	{

		fprintf(finit,"%d\n",npart);
		fflush(finit);

		fprintf(finit,"teste\n");
		fflush(finit);

		for(i=0;i<l3;++i)
	    {
			//vizinhos
			site = i;
			x = site%l;
			y = (site/l)%l;
			z = site/l2;

			if ( (s[i] == 1) ) //Água
			{
				fprintf(finit,"O  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do IF
	
			else if ( (s[i] == 2) ) //óleo
			{
				fprintf(finit,"V  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do ELSE IF

			else if ( (s[i] == 9 ) ) //óleo
			{
				fprintf(finit,"C  %d  %d  %d\n",x,y,z);
				fflush(finit);
			}//fim do ELSE IF

		}//Fim do for

	}//fim do IF

	else if ( (m == 1) )
	{

		npart = npart;
		fprintf(fxyz,"%d\n",npart);
		fflush(fxyz);

		fprintf(fxyz,"teste\n");
		fflush(fxyz);

		for(i=0;i<l3;++i)
	    {
			//vizinhos
			site = i;
			x = site%l;
			y = (site/l)%l;
			z = site/l2;

			if ( (s[i] == 1) ) //Água
			{
				fprintf(fxyz,"O  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do IF
	
			else if ( (s[i] == 2) ) //óleo
			{
				fprintf(fxyz,"V  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do ELSE IF

			else if ( (s[i] == 9) ) //óleo
			{
				fprintf(fxyz,"C  %d  %d  %d\n",x,y,z);
				fflush(fxyz);
			}//fim do ELSE IF

		}//Fim do for

	}//fim do ELSE IF

 return;

}

/*********************************************************************************************
*                               Rotina que faz a dinâmica do MC                              *
*********************************************************************************************/
void dynamics(int *s,int num_steps, double Aw)
{
	int     i, j, k, hs, xi, xf, dx, yi, yf, dy;
/*	int 	z; */
	int     site,s_site, s_teste, label,soma; 
	int     nw,ns,no,ng;
	double	x_CM, y_CM, z_CM;
	double  temp_e, delta_s, delta_v, delta_g, delta_o, delta_m, delta;
	long int count=0, count_w=0;
	int x_just_air = 0, y_just_air = 0, x_pos, y_pos;

	x_CM = 0.0;
	y_CM = 0.0;
	z_CM = 0.0;
	// Iterate over yz planes to find the first x plane without water
	for(i=0; i<l; i++) 
	{
		x_just_air = i;
		for(j=0;j<l;j++) 
		{
			for(k=h_base; k<l; k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					x_just_air = -1;
					break; // Exit the inner loop
				}
			}
			if(x_just_air == -1) break; // Exit the outer loop
		}
		if(x_just_air > -1) break; // Exit the outermost loop
	}

	// Iterate over xz planes to find the first y plane without water
	for(j=0; j<l; j++) 
	{
		y_just_air = j;
		for(i=0;i<l;i++) 
		{
			for(k=h_base; k<l; k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					y_just_air = -1;
					break; // Exit the inner loop
				}
			}
			if(y_just_air == -1) break; // Exit the outer loop
		}
		if(y_just_air > -1) break; // Exit the outermost loop
	}

	for(i=0;i<=l;i++) 
	{
		for(j=0;j<=l;j++) 
		{
			for(k=h_base;k<l;k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					count_w++;

					if(i < x_just_air){
						x_pos = i + l;
					} else
					{
						x_pos = i;
					}
					if(j < y_just_air){
						y_pos = j + l;
					} else
					{
						y_pos = j;
					}
					x_CM = x_CM + x_pos;
					y_CM = y_CM + y_pos;
					z_CM = z_CM + k;
            	}
			}
		}
	}
	// Calculate center of mass and shift back to original range
	x_CM = fmod((float) x_CM/(float) count_w,l);
	y_CM = fmod((float) y_CM/(float) count_w,l);
	z_CM = fmod((float) z_CM/(float) count_w,l);

	for(j = 0; j < t_vol ; ++j)
    {

//--------------------------------------------------------------------------------------------
//   Sorteando um sítio na interface
//--------------------------------------------------------------------------------------------	

		label = (int) (FRANDOM*int_label);
		site = w_inter[label];
		hs = (site/l2-h_base);
		
/*		z = site/l2;*/
		s_site=s[site];

		trial++;
      
		acc_interface++;
      
		acc_teste++;

		delta_s=0;
		delta_v=0;
		delta_o=0;
		delta_g=0;
		delta_m=0;
      
		nw=nv_w[site];
		no=nv_o[site];
		ns=nv_s[site];
		ng=nv_g[site];

//--------------------------------------------------------------------------------------------
//   Determinando qual será a troca
//--------------------------------------------------------------------------------------------	

		if(f_o == 0)
		{

			if(s_site==0)
			{
				s_teste = 1; //troca de ar-agua
			}//fim do IF
			else
			{
				s_teste = 0; //troca de agua-ar
			}//fim do ELSE

		}//fim do IF
		else if(f_o == 1)
		{

			if(s_site==0)
			{
				s_teste = 2; //troca de ar-oleo
			}//fim do IF
			else
			{
				s_teste = 0; //troca de oleo-ar
			}//fim do ELSE

		}//fim do ELSE IF
		else
		{
			s_teste = s_site;
			while (s_teste == s_site)
			{
				s_teste = (int) (FRANDOM*3);
			}//fim do WHILE
			
			if(s_teste==s_site) printf("ERRO: s_teste=s_site!!\n");

		}//fim do ELSE

		soma = s_site +s_teste;	

//--------------------------------------------------------------------------------------------
//   Calculando a energia da troca
//--------------------------------------------------------------------------------------------
		
		if (soma==1)  //água e ar
		{

			if (s_site==0)
			{
				xf = site % l;
				dx = xf - x_CM;
				// Adjust dx for periodic boundary conditions
				if (dx > l/2) {
					dx -= l;
				} else if (dx < -l/2) {
					dx += l;
				}
				xi = (dx > 0) ? (xf - 1) : (xf + 1);
				dx = xf - xi;

				yf = (site/l) % l;
				dy = yf - y_CM;
				// Adjust dx for periodic boundary conditions
				if (dy > l/2) {
					dy -= l;
				} else if (dy < -l/2) {
					dy += l;
				}
				yi = (dy > 0) ? (yf - 1) : (yf + 1);
				dy = yf - yi;
				
				// Set delta_m based on the sign of dx
				delta_m = -Aw * (1*dx + 1*dy); // Displacement dot versor (x=1,y=1), is this rigth?
				delta_g = Gw*hs;
				delta_s = (ng-nw)*eps_WG + ns*(eps_SW-eps_SG) + (eps_WO-eps_OG)*no;
/*				delta_v = (1+2*(vol-t_vol))*Lambda_w; // Ganho um  líquido*/
				delta_v = (1+2*(vol_w-t_vol_w))*Lambda_w; // Ganha uma água

			}//fim do IF

			else //if (s_site==1)
			{
				xi = site % l;
				dx = x_CM - xi;
				// Adjust dx for periodic boundary conditions
				if (dx > l/2) {
					dx -= l;
				} else if (dx < -l/2) {
					dx += l;
				}
				xf = (dx > 0) ? (xi + 1) : (xi - 1);
				dx = xf - xi;

				yi = (site/l) % l;
				dy = y_CM - yi;
				// Adjust dx for periodic boundary conditions
				if (dy > l/2) {
					dy -= l;
				} else if (dy < -l/2) {
					dy += l;
				}
				yf = (dy > 0) ? (yi + 1) : (yi - 1);
				dy = yf - yi;
				
				// Set delta_m based on the sign of dx
				delta_m = -Aw * (1*dx + 1*dy); // Displacement vector dot versor (x=1,y=1), is this rigth?
				delta_g = -Gw*hs; 
				delta_s = (nw-ng)*eps_WG + ns*(eps_SG-eps_SW) + (eps_OG-eps_WO)*no;
/*				delta_v = (1-2*(vol-t_vol))*Lambda_w; // Perco um líquido*/
				delta_v = (1-2*(vol_w-t_vol_w))*Lambda_w; // Perco uma água

			}// fim do ELSE
			
		} //fim do IF

		if ((soma==2) && (s_teste!=s_site))  //óleo e ar
		{

			if (s_site==0)
			{

				delta_g = Go*hs;
				delta_s = (ng-no)*eps_OG + ns*(eps_SO-eps_SG) + (eps_WO-eps_WG)*nw;
/*				delta_v = (1+2*(vol-t_vol))*Lambda_w; //Ganho um líquido */
				delta_o = (1+2*(vol_o-t_vol_o))*(Lambda_o); // Ganho um um óleo
	  

			}//fim do IF

			else //if (s_site==2)
			{
	
				delta_g = -Go*hs; 
				delta_s = (no-ng)*eps_OG + ns*(eps_SG-eps_SO) + (eps_WG-eps_WO)*nw;
/*				delta_v = (1-2*(vol-t_vol))*Lambda_w; // Perco um líquido */
				delta_o = (1-2*(vol_o-t_vol_o))*Lambda_o;  // Perco um  um óleo

						  

			}// fim do ELSE

		} //fim do  IF

		if (soma==3)  //óleo e água
		{
			
			if (s_site==1)
			{	
				delta_g = (Go-Gw)*hs;
				delta_s = (nw-no)*eps_WO + (eps_OG-eps_WG)*ng + (eps_SO-eps_SW)*ns;
				delta_o = (1+2*(vol_o-t_vol_o))*Lambda_o; //Ganho um óleo
				delta_v = (1-2*(vol_w-t_vol_w))*Lambda_w; // Perco uma água
				

			}//fim do IF

			else //if (s_site==2)
			{
				delta_g = (Gw-Go)*hs; 
				delta_s = (no-nw)*eps_WO + (eps_WG-eps_OG)*ng + (eps_SW-eps_SO)*ns;
				delta_o = (1-2*(vol_o-t_vol_o))*Lambda_o;  // Perco um óleo
				delta_v = (1+2*(vol_w-t_vol_w))*Lambda_w; // Ganha uma água

			}// fim do ELSE	

		}//fim do IF

		if (soma>3)
		{
			printf("ENTROU!!\n");
			delta_g = 100000.0; 
			delta_s = 100000.0;
			delta_v = 100000.0;
		}//fim do else
 
		delta=delta_v+delta_s+delta_g+delta_o+delta_m;

//--------------------------------------------------------------------------------------------
//  Testando a troca
//--------------------------------------------------------------------------------------------
      
		if (delta<=0) 
		{

/*			if(num_steps > n_teste)*/
/*			{*/
/*				if( z > (h+h_base+10) && s_teste==2 ) teste_flip_droplet(site,s_teste,delta_s, delta);*/
/*				if( z < (h+h_base+1) && s_site==2 ) teste_flip_surface(site,s_teste,delta_s, delta);*/
/*			}*/

			flip_spin(site,s_teste);
			count++;
			accepted++;

		}//fim do IF

 		else 
		{
			if (T > 0.0) 
			{
				temp_e=FRANDOM;
				if(temp_e<exp(-delta/(kB*T))) 
				{

/*					if(num_steps > n_teste)*/
/*					{*/
/*						if( z > (h+h_base+10) && s_teste==2 ) teste_flip_droplet(site,s_teste,delta_s, delta);*/
/*						if( z < (h+h_base+1) && s_site==2 ) teste_flip_surface(site,s_teste,delta_s, delta);*/
/*					}*/

					flip_spin(site,s_teste);
					count++;
					accepted++;

				}//fim do IF
			}//fim do IF
		}//dim do ELSE

	} //fim do FOR

  return;
}
/*********************************************************************************************
*                               Rotina que flipa o spin do sítio                             *
*********************************************************************************************/
void flip_spin (long site, long new_s)
{
	int i,neighbour;
	int soma;
	int neigh[26];
	int x,y,z,xn,yn,zn;

//--------------------------------------------------------------------------------------------
//  Atualizando o volume
//--------------------------------------------------------------------------------------------

	soma = s[site]+new_s;

	if (soma==1)  //água e ar
	{
		if (s[site]==0)
		{
			++vol;
			++vol_w;
		}//fim do IF

		else //if (s_site==1)
		{
			--vol;
			--vol_w;
		}// fim do ELSE		
	} //fim do IF

	else if (soma==2)  //óleo e ar
	{
		if (s[site]==0)
		{
			++vol;
			++vol_o;
		}//fim do IF

		else //if (s_site==1)
		{
			--vol;
			--vol_o;
		}// fim do ELSE		

	} //fim do ELSE IF

	else if (soma==3)  //óleo e água
	{
		if (s[site]==1)
		{
			++vol_o;
			--vol_w;
		}//fim do IF

		else //if (s_site==2)
		{
			--vol_o;
			++vol_w;
		}// fim do ELSE					
					
	}//fim do ELSE IF

	n_w = (int) vol_w;
	n_o = (int) vol_o;

	
//--------------------------------------------------------------------------------------------
//  Atualizando o numero de vizinhos
//--------------------------------------------------------------------------------------------
  
	x = site%l;
	y = (site/l)%l;
	z = site/l2;
      
	// 0
	yn = (y==0)?l-1:y-1;
	neigh[0]=x + yn*l + z*l2;
	// 1
	yn = (y==l-1)?0:y+1;
	neigh[1] = x + yn*l + z*l2;
	// 2
	xn = (x==l-1)?0:x+1;
	neigh[2] = xn + y*l + z*l2;
	// 3
	xn = (x==0)?l-1:x-1;
	neigh[3] = xn + y*l + z*l2;
	// 4
	zn = (z==0)?l-1:z-1;
	neigh[4] = x + y*l + zn*l2;
	// 5
	zn = (z==l-1)?0:z+1;
	neigh[5] = x + y*l + zn*l2;
	//6 
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	neigh[6] = xn + yn*l + z*l2;
	// 7
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	neigh[7] = xn + yn*l + z*l2;
	// 8 
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	neigh[8] = xn + yn*l + z*l2;
	// 9
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	neigh[9] = xn + yn*l + z*l2;
	// 10
	xn = (x==l-1)?0:x+1;
	zn = (z==0)?l-1:z-1;
	neigh[10] = xn + y*l + zn*l2;
	// 11
	xn = (x==0)?l-1:x-1;
	zn = (z==0)?l-1:z-1;
	neigh[11] = xn + y*l + zn*l2;
	// 12
	xn = (x==l-1)?0:x+1;
	zn = (z==l-1)?0:z+1;
	neigh[12] = xn + y*l + zn*l2;
	// 13
	xn = (x==0)?l-1:x-1;
	zn = (z==l-1)?0:z+1;
	neigh[13] = xn + y*l + zn*l2;
	// 14
	zn = (z==0)?l-1:z-1;
	yn = (y==0)?l-1:y-1;
	neigh[14] = x + yn*l + zn*l2;
	// 15
	zn = (z==l-1)?0:z+1;
	yn = (y==0)?l-1:y-1;
	neigh[15] = x + yn*l + zn*l2;
	// 16
	zn = (z==0)?l-1:z-1;
	yn = (y==l-1)?0:y+1;
	neigh[16] = x + yn*l + zn*l2;
	// 17
	zn = (z==l-1)?0:z+1;
	yn = (y==l-1)?0:y+1;
	neigh[17] = x + yn*l + zn*l2;
	// 18
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	zn = (z==0)?l-1:z-1;
	neigh[18] = xn + yn*l + zn*l2;
	// 19
	xn = (x==l-1)?0:x+1;
	yn = (y==0)?l-1:y-1;
	zn = (z==l-1)?0:z+1;
	neigh[19] = xn + yn*l + zn*l2;
	// 20
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	zn = (z==0)?l-1:z-1;
	neigh[20] = xn + yn*l + zn*l2;
	// 21
	xn = (x==0)?l-1:x-1;
	yn = (y==0)?l-1:y-1;
	zn = (z==l-1)?0:z+1;
	neigh[21] = xn + yn*l + zn*l2;
	// 22
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	zn = (z==0)?l-1:z-1;
	neigh[22] = xn + yn*l + zn*l2;
	// 23
	xn = (x==l-1)?0:x+1;
	yn = (y==l-1)?0:y+1;
	zn = (z==l-1)?0:z+1;
	neigh[23] = xn + yn*l + zn*l2;
	// 24
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	zn = (z==0)?l-1:z-1;
	neigh[24] = xn + yn*l + zn*l2;
	// 25	  
	xn = (x==0)?l-1:x-1;
	yn = (y==l-1)?0:y+1;
	zn = (z==l-1)?0:z+1;
	neigh[25] = xn + yn*l + zn*l2;

	for(i=0;i<t_close_neigh;++i)
    {
		neighbour=neigh[i];

		if (soma==1)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				++nv_cw[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				--nv_cw[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do IF

		else if (soma==2)
		{
			if (new_s==2)
			{
				++nv_o[neighbour];
				++nv_co[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_o[neighbour];
				--nv_co[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

		else if (soma==3)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				++nv_cw[neighbour];
				--nv_o[neighbour];
				--nv_co[neighbour];

			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				--nv_cw[neighbour];
				++nv_o[neighbour];
				++nv_co[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

	}//fim do FOR


	for (i=t_close_neigh;i<t_neigh ;++i)
    {
		neighbour=neigh[i];

		if (soma==1)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do IF

		else if (soma==2)
		{
			if (new_s==2)
			{
				++nv_o[neighbour];
				--nv_g[neighbour];
			} //fim do IF
			else 
			{
				--nv_o[neighbour];
				++nv_g[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

		else if (soma==3)
		{
			if (new_s==1)
			{
				++nv_w[neighbour];
				--nv_o[neighbour];

			} //fim do IF
			else 
			{
				--nv_w[neighbour];
				++nv_o[neighbour];
			}//fim do ELSE
		}//fim do ELSE IF

    }//fim do FOR

//--------------------------------------------------------------------------------------------
//  Atualizando sítios na interface
//--------------------------------------------------------------------------------------------

	for (i=0;i<t_close_neigh;++i)
	{
		neighbour=neigh[i];

		if (new_s==0)
		{
			if ( (s[neighbour]==1) || (s[neighbour]==2) )
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}

	    	}//fim do IF

			else if ((s[neighbour]==0))
			{
	   			if ((nv_cw[neighbour]==0) && (nv_co[neighbour]==0) && (inter_pos[neighbour]==1))
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF

		}//fim do IF
		
		else if(new_s==1)
		{
			if (s[neighbour]==0 || (s[neighbour]==2))
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}
			}//fim do IF

			else if (s[neighbour]==1)
			{

			    if ((nv_cw[neighbour]==t_close_neigh) && (inter_pos[neighbour]==1) )
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF
		}//fim do ELSE IF

		else if(new_s==2)
		{
			if (s[neighbour]==0 || (s[neighbour]==1))
			{
				if (inter_pos[neighbour]==0) {add_to_interface(neighbour);}
				if (inter_pos[site]==0) {add_to_interface(site);}
			}//fim do IF

			else if (s[neighbour]==2)
			{

			    if ((nv_co[neighbour]==t_close_neigh) && (inter_pos[neighbour]==1) )
				{
					remove_from_interface(neighbour);
				} //fim do IF
			}//fim do ELSE IF
		}//fim do ELSE IF

	}//fim do FOR
      
	if (new_s==1) 
	{
		if ((nv_cw[site]==t_close_neigh)&&(inter_pos[site]==1)) {remove_from_interface(site);}
    }//fim do IF
	else if (new_s==2)
	{
		if ((nv_co[site]==t_close_neigh)&&(inter_pos[site]==1)) {remove_from_interface(site);}
	}//fim do ELSE IF
  	else if (new_s==0)
    {
      if (((nv_co[site]+nv_cw[site])==0)&&(inter_pos[site]==1)) {remove_from_interface(site);}
    }//fim do ELSE IF
  
	s[site]=new_s;
     
  return;
}
/*********************************************************************************************
*                      Rotina que adiciona um sítio da lista da interface                    *
*********************************************************************************************/
void add_to_interface(int site)
{
 
	inter_pos[site]=1;
	inter_mtx[site]=int_label;
	w_inter[int_label]=site;
	++int_label;

  return;
}
/*********************************************************************************************
*                       Rotina que remove um sítio da lista da interface                     *
*********************************************************************************************/
void remove_from_interface(int site)
{
	int site_label, last_interface_site;
  
	inter_pos[site]=0;

  	site_label=inter_mtx[site];
  	inter_mtx[site]=-1;
  
	last_interface_site=w_inter[int_label-1];
  	w_inter[site_label]=last_interface_site;
	inter_mtx[last_interface_site]=site_label;
  
	w_inter[int_label-1]=-1;
	--int_label;
  
  return;
}

/*********************************************************************************************
*                               Rotina que abre os ponteiros                                 *
*********************************************************************************************/
void create_pointers(void)
{
  s=create_int_pointer(l3);
  nv_cw=create_int_pointer(l3);
  nv_co=create_int_pointer(l3);

  nv_w=create_int_pointer(l3);
  nv_o=create_int_pointer(l3);
  nv_s=create_int_pointer(l3);
  nv_g=create_int_pointer(l3);

  inter_pos=create_int_pointer(l3);
  w_inter=create_int_pointer(l3);
  inter_mtx=create_int_pointer(l3);

  return;
}


/*********************************************************************************************
*                                Rotina que libera os ponteiros                              *
*********************************************************************************************/
void free_pointers(void)
{
	free(s);
	free(nv_cw);
	free(nv_co);
	free(nv_w);
	free(nv_o);
	free(nv_s);
	free(nv_g);
	free(inter_pos);
  	free(w_inter);
	free(inter_mtx);

  return;
}

/*********************************************************************************************
*                                Rotina que calcula observáveis                              *
*********************************************************************************************/
void measure_angle(int num_steps) 
{
  
	int i,j,k;
	int xmin,xmax,ymin,ymax;
	int xmin_o,xmax_o,ymin_o,ymax_o;
	int h_film, h_film_o;
	double H2_B2_x,H2_B2_y;
	double H2_B2_x_o,H2_B2_y_o;
	int pil_xmin,pil_xmax,pil_ymin,pil_ymax;
	int num_pil_x,num_pil_y;
	int num_pil_x_o,num_pil_y_o;
	int height,site;
	double volume_below, volume_below_o;
	double volume_res, volume_res_o;
	int vbt, vrt;
	double rx,ry;
	double rx_o,ry_o;
	int count_w, count_o;
	int count_nw, count_no;
	int nw, no;
	int int_w, int_o;
	double per_w, per_o;
	double frac_w, frac_o;
	int htest;
	double vol_p_w,vol_p_o;
	double x_CM, y_CM, z_CM;
	int x_just_air , y_just_air , x_pos, y_pos;

// -------------------------------------------------------------------------------------------
// Determinando a altura do filme para água e óleo

	h_film = 0;
	h_film_o = 0;
	count_w=0;
	count_o=0;

	for(i=0;i<20;++i)
		for(j=0;j<20;++j)
			for(k=h+h_base;k<l;++k)
	{

				site = i +j*l + k*l2;
				nw=nv_w[site];
				no=nv_o[site];

				// para água
				if ((s[site]==1) && (nw>2)) 
				{
					++count_w;
					if (k>h_film) h_film=k;
				} // fim do IF

				// para óleo
				if ((s[site]==2) && (no>2)) 
				{
					++count_o;
					if (k>h_film_o) h_film_o=k;
				} // fim do IF			
		
	} //fim do FOR


  	h_film-= h+h_base;
	h_film_o-= h+h_base;

	if (count_w==0)  h_film=0;
	if (count_o==0)  h_film_o=0;

// -------------------------------------------------------------------------------------------
// Determinando quantos sítios da interface são agua e óleo

	int_w = 0;
	int_o = 0;

	for(i=0;i<int_label;++i)
    {
		site = w_inter[i];

		if (s[site]==1) ++int_w;
		if (s[site]==2) ++int_o;

	}

// -------------------------------------------------------------------------------------------
// Determinando a altura da gota para água e óleo
// Nesse calculo só entram os sítios que não participam do FILME

	count_w=0;
	count_o=0;

	drop_height=0;
	drop_bottom=l;
	drop_height_o=0;
	drop_bottom_o=l;

	for(i=0;i<int_label;++i)
    {
		site = w_inter[i];
		height = site/l2;
		nw=nv_w[site];
		no=nv_o[site];
		
		// Só pega sítios de água que não estão no reservatório NEM NOS POROS
		if ( (s[site]==1) && (height>h+h_base) && (nw>=4)) 
		{

			++count_w;
			if (height>drop_height) drop_height=height;
			if (height<drop_bottom) drop_bottom=height;

		}//fim do IF

		// Só pega sítios de óleo que não estão no reservatório  NEM NOS POROS	
		if ( (s[site]==2) && (height>h+h_base) && (no>=4)) 
		{

			++count_o;
			if (height>drop_height_o) drop_height_o=height;
			if (height<drop_bottom_o) drop_bottom_o=height;

		}//fim do IF


    } //fim do FOR

  	drop_height-= h+h_base;
	drop_height_o-= h+h_base;

	if (count_w==0)  drop_height=0;
	if (count_o==0)  drop_height_o=0;	

	per_w = (float)count_w/int_w;
	per_o = (float)count_o/int_o;	

// -------------------------------------------------------------------------------------------
// Determinando numero de vizinhos de mesmo tipo em z=h+hbase
	
	count_w=0;
	count_o=0;
	count_nw=0;
	count_no=0;		

	for(i=0;i<l;++i)
		for(j=0;j<l;++j)
	{

		site = (h+h_base)*l2+j*l+i;
		nw=nv_w[site];
		no=nv_o[site];
		
		if (s[site]==1) 
		{
			++count_w;
			count_nw += nw;		
		}

		if (s[site]==2) 
		{
			++count_o;
			count_no += no;		
		}

	}//fim do FOR

	frac_w = (float)count_nw/count_w;
	frac_o = (float)count_no/count_o;

// -------------------------------------------------------------------------------------------
// Determinando o x e y minimo e maximo da gota para água e óleo

	major_axis_x=0,major_axis_y=0;
	major_axis_x_o=0,major_axis_y_o=0;

	xmin=l;
	xmax=0;
	ymin=l;
	ymax=0;

	xmin_o=l;
	xmax_o=0;
	ymin_o=l;
	ymax_o=0;

	count_w=0;
	count_o=0;

	for(i=0;i<l;++i)
		for(j=0;j<l;++j)
	{

		site = (h+h_base)*l2+j*l+i;
		nw=nv_w[site];
		no=nv_o[site];

		if ((s[site]==1) && (nw>=2))
		{
			++count_w;
			if (i<xmin) xmin=i;
			if (i>xmax) xmax=i;
			if (j<ymin) ymin=j;
			if (j>ymax) ymax=j;

		} //fim do IF

		if (s[site]==2 && (no>=2))
		{
			++count_o;
			if (i<xmin_o) xmin_o=i;
			if (i>xmax_o) xmax_o=i;
			if (j<ymin_o) ymin_o=j;
			if (j>ymax_o) ymax_o=j;
		} //fim do IF

	} //fim dos FOR

	if (count_w==0)  
	{
		xmin=0;
		xmax=0;
		ymin=0;
		ymax=0;
	} // fim do IF

	if (count_o==0)  
	{
		xmin_o=0;
		xmax_o=0;
		ymin_o=0;
		ymax_o=0;
	} // fim do IF

// -------------------------------------------------------------------------------------------
// Determinando número de pilares que a gota toca - água e óleo

	pil_xmin = (int)xmin/(w+a);
	pil_xmax = (int)xmax/(w+a);
	pil_ymin = (int)ymin/(w+a);
	pil_ymax = (int)ymax/(w+a);

	num_pil_x = pil_xmax - pil_xmin;
	num_pil_y = pil_ymax - pil_ymin;

	pil_xmin = (int)xmin_o/(w+a);
	pil_xmax = (int)xmax_o/(w+a);
	pil_ymin = (int)ymin_o/(w+a);
	pil_ymax = (int)ymax_o/(w+a);

	num_pil_x_o = pil_xmax - pil_xmin;
	num_pil_y_o = pil_ymax - pil_ymin;

// -------------------------------------------------------------------------------------------
// Determinando o raio da base da gota em x e y - água e óleo

	base_radius_x = 0.5*(xmax-xmin);
	base_radius_y = 0.5*(ymax-ymin);
	centre_x =(int) (0.5*(xmax+xmin));
	centre_y =(int) (0.5*(ymax+ymin));

	base_radius_x_o = 0.5*(xmax_o-xmin_o);
	base_radius_y_o = 0.5*(ymax_o-ymin_o);
	centre_x_o =(int) (0.5*(xmax_o+xmin_o));
	centre_y_o =(int) (0.5*(ymax_o+ymin_o));
// -------------------------------------------------------------------------------------------
// Caclulando as hipotenusas - água e óleo

	H2_B2_x=drop_height*drop_height+base_radius_x*base_radius_x;
	H2_B2_y=drop_height*drop_height+base_radius_y*base_radius_y;

	H2_B2_x_o=drop_height_o*drop_height_o+base_radius_x_o*base_radius_x_o;
	H2_B2_y_o=drop_height_o*drop_height_o+base_radius_y_o*base_radius_y_o;

// -------------------------------------------------------------------------------------------
// Raio da gota em x e y - água e óleo

	rx = H2_B2_x/(2.0*drop_height);
	if(rx>(l/2)) rx=(float)l/2;

	ry = H2_B2_y/(2.0*drop_height);
	if(ry>(l/2)) ry=(float)l/2;

	drop_radius_x = rx;
	drop_radius_y = ry;

	rx_o = H2_B2_x_o/(2.0*drop_height_o);
	if(rx_o>(l/2)) rx_o=(float)l/2;

	ry_o = H2_B2_y_o/(2.0*drop_height_o);
	if(ry_o>(l/2)) ry_o=(float)l/2;

	// preciso pra visualizacao
	drop_radius_x_o = rx_o;
	drop_radius_y_o = ry_o;

// ------------------------------------------------------------------------------------------- 
// angulos - água  e óleo

	angle_x = asin(base_radius_x/rx);
	angle_y = asin(base_radius_y/ry);


	if (drop_height>rx) 
	{
		angle_x = M_PI - angle_x;
	} // fim do IF

	if (drop_height>ry) 
	{
		angle_y = M_PI - angle_y;
	}// fim do IF

	angle_x = angle_x*180.0/M_PI; 
	angle_y = angle_y*180.0/M_PI;
 

	angle_x_o = asin(base_radius_x_o/rx_o);
	angle_y_o = asin(base_radius_y_o/ry_o);

	if (drop_height_o>rx_o) 
	{
		angle_x_o = M_PI - angle_x_o;
	} // fim do IF

	if (drop_height_o>ry_o) 
	{
		angle_y_o = M_PI - angle_y_o;
	}// fim do IF

	angle_x_o = angle_x_o*180.0/M_PI ;
	angle_y_o = angle_y_o*180.0/M_PI; 

	if(num_steps > (mc_steps/2))
	{
		thmed = thmed+(0.5*angle_x+0.5*angle_y);
		++c_theta;
	}//fim do IF

// ------------------------------------------------------------------------------------------- 
// volume entre os pilares - água e óleo

	volume_below=0;
	volume_below_o=0;
	vbt=0;

	for(i=0;i<=l;i++) 
	{
		for(j=0;j<=l;j++) 
		{
      		for(k=h_base;k<(h_base+h);k++) 
			{

				vbt++;
				site = i +j*l + k*l2;

				if( (s[site]==1) ) volume_below++;
				if( (s[site]==2) ) volume_below_o++;

			}//fim do FOR
		}//fim do FOR
	}//fim do FOR

	volume_below = volume_below/n_w;
	volume_below_o = volume_below_o/n_o;

// ------------------------------------------------------------------------------------------- 
// volume no reservatório - água e óleo

	volume_res=0;
	volume_res_o=0;
	vrt=0;

	for(i=0;i<=l;i++) 
	{
		for(j=0;j<=l;j++) 
		{
      		for(k=0;k<(h_base);k++) 
			{
			
				vrt++;
				site = i +j*l + k*l2;

				if( (s[site]==1) ) volume_res++;
				if( (s[site]==2) ) volume_res_o++;

			}//fim do FOR
		}//fim do FOR
	}//fim do FOR


	volume_res = volume_res/n_w;
	volume_res_o = volume_res_o/n_o;
// ------------------------------------------------------------------------------------------- 
// Porcentagem de água e óleo na gota

	count_w=0;
	count_o=0;

	if(h_film==0 && h_film_o==0) 
	{
		htest=h+h_base;

	} //fim do IF
	else if (h_film <= h_film_o)
	{
		htest=h+h_base + h_film_o;
	}//fim do ELSE IF
	else
	{
		htest=h+h_base + h_film;
	} //fim do IF

	for(i=0;i<=l;i++) 
	{
		for(j=0;j<=l;j++) 
		{
      		for(k=htest;k<l;k++) 
			{
				site = i +j*l + k*l2;

				if( (s[site]==1) ) count_w++;
				if( (s[site]==2) ) count_o++;

			}//fim do FOR
		}//fim do FOR
	}//fim do FOR

	vol_p_w = (float)count_w/(count_w+count_o);
	vol_p_o = (float)count_o/(count_w+count_o);


	if(num_steps > (mc_steps/2))
	{
		fvmed_w = fvmed_w+vol_p_w;
		++c_fv_w;

		fvmed_o = fvmed_o+vol_p_o;
		++c_fv_o;
	}//fim do IF

// -------------------------------------------------------------------------------------------
// Rotina que calcula o centro de massa para uma gota de água
	x_CM = 0.0;
	y_CM = 0.0;
	z_CM = 0.0;
	count_w = 0;
	// Iterate over yz planes to find the first x plane without water
	for(i=0; i<l; i++) 
	{
		x_just_air = i;
		for(j=0;j<l;j++) 
		{
			for(k=h_base; k<l; k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					x_just_air = -1;
					break; // Exit the inner loop
				}
			}
			if(x_just_air == -1) break; // Exit the outer loop
		}
		if(x_just_air > -1) break; // Exit the outermost loop
	}

	// Iterate over xz planes to find the first y plane without water
	for(j=0; j<l; j++) 
	{
		y_just_air = j;
		for(i=0;i<l;i++) 
		{
			for(k=h_base; k<l; k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					y_just_air = -1;
					break; // Exit the inner loop
				}
			}
			if(y_just_air == -1) break; // Exit the outer loop
		}
		if(y_just_air > -1) break; // Exit the outermost loop
	}

	for(i=0;i<=l;i++) 
	{
		for(j=0;j<=l;j++) 
		{
			for(k=h_base;k<l;k++) 
			{
				site = i +j*l + k*l2;
				if( (s[site]==1) ) 
				{
					count_w++;

					if(i < x_just_air){
						x_pos = i + l;
					} else
					{
						x_pos = i;
					}
					if(j < y_just_air){
						y_pos = j + l;
					} else
					{
						y_pos = j;
					}
					x_CM = x_CM + x_pos;
					y_CM = y_CM + y_pos;
					z_CM = z_CM + k;
            	}
			}
		}
	}
	// Calculate center of mass and shift back to original range
	x_CM = fmod((float) x_CM/(float) count_w,l);
	y_CM = fmod((float) y_CM/(float) count_w,l);
	z_CM = fmod((float) z_CM/(float) count_w,l);
// -------------------------------------------------------------------------------------------
  
	fprintf(fp1,"%8d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", num_steps, (int)vol, (int)vol_w, (int)vol_o, calculate_energy(), base_radius_x, base_radius_y, base_radius_x_o, base_radius_y_o, rx, ry, rx_o, ry_o, angle_x, angle_y, angle_x_o, angle_y_o, num_pil_x,num_pil_y, num_pil_x_o,num_pil_y_o, volume_below, volume_below_o, volume_res,volume_res_o, per_w, per_o, frac_w, frac_o, vol_p_w, vol_p_o, x_CM, y_CM, z_CM);
	fflush(fp1);

  return;
}


/*********************************************************************************************
*                                 Rotina que calcula a energia                               *
*********************************************************************************************/
double calculate_energy(void)
{
	double total=0.0;
	int i,j;
	int neigh[26];  
	int x,y,z,xn,yn,zn,site;
	long int num_interface=0;
	double ene_interface=0.0,ene_grav=0.0;
   
	for(i=0;i<l3;++i)
    {
		//vizinhos
		site = i;
		x = site%l;
		y = (site/l)%l;
		z = site/l2;
      
		// 0
		yn = (y==0)?l-1:y-1;
		neigh[0]=x + yn*l + z*l2;
		// 1
		yn = (y==l-1)?0:y+1;
		neigh[1] = x + yn*l + z*l2;
		// 2
		xn = (x==l-1)?0:x+1;
		neigh[2] = xn + y*l + z*l2;
		// 3
		xn = (x==0)?l-1:x-1;
		neigh[3] = xn + y*l + z*l2;
		// 4
		zn = (z==0)?l-1:z-1;
		neigh[4] = x + y*l + zn*l2;
		// 5
		zn = (z==l-1)?0:z+1;
		neigh[5] = x + y*l + zn*l2;
		//6 
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		neigh[6] = xn + yn*l + z*l2;
		// 7
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		neigh[7] = xn + yn*l + z*l2;
		// 8 
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		neigh[8] = xn + yn*l + z*l2;
		// 9
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		neigh[9] = xn + yn*l + z*l2;
		// 10
		xn = (x==l-1)?0:x+1;
		zn = (z==0)?l-1:z-1;
		neigh[10] = xn + y*l + zn*l2;
		// 11
		xn = (x==0)?l-1:x-1;
		zn = (z==0)?l-1:z-1;
		neigh[11] = xn + y*l + zn*l2;
		// 12
		xn = (x==l-1)?0:x+1;
		zn = (z==l-1)?0:z+1;
		neigh[12] = xn + y*l + zn*l2;
		// 13
		xn = (x==0)?l-1:x-1;
		zn = (z==l-1)?0:z+1;
		neigh[13] = xn + y*l + zn*l2;
		// 14
		zn = (z==0)?l-1:z-1;
		yn = (y==0)?l-1:y-1;
		neigh[14] = x + yn*l + zn*l2;
		// 15
		zn = (z==l-1)?0:z+1;
		yn = (y==0)?l-1:y-1;
		neigh[15] = x + yn*l + zn*l2;
		// 16
		zn = (z==0)?l-1:z-1;
		yn = (y==l-1)?0:y+1;
		neigh[16] = x + yn*l + zn*l2;
		// 17
		zn = (z==l-1)?0:z+1;
		yn = (y==l-1)?0:y+1;
		neigh[17] = x + yn*l + zn*l2;
		// 18
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		zn = (z==0)?l-1:z-1;
		neigh[18] = xn + yn*l + zn*l2;
		// 19
		xn = (x==l-1)?0:x+1;
		yn = (y==0)?l-1:y-1;
		zn = (z==l-1)?0:z+1;
		neigh[19] = xn + yn*l + zn*l2;
		// 20
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		zn = (z==0)?l-1:z-1;
		neigh[20] = xn + yn*l + zn*l2;
		// 21
		xn = (x==0)?l-1:x-1;
		yn = (y==0)?l-1:y-1;
		zn = (z==l-1)?0:z+1;
		neigh[21] = xn + yn*l + zn*l2;
		// 22
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		zn = (z==0)?l-1:z-1;
		neigh[22] = xn + yn*l + zn*l2;
		// 23
		xn = (x==l-1)?0:x+1;
		yn = (y==l-1)?0:y+1;
		zn = (z==l-1)?0:z+1;
		neigh[23] = xn + yn*l + zn*l2;
		// 24
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		zn = (z==0)?l-1:z-1;
		neigh[24] = xn + yn*l + zn*l2;
		// 25	  
		xn = (x==0)?l-1:x-1;
		yn = (y==l-1)?0:y+1;
		zn = (z==l-1)?0:z+1;
		neigh[25] = xn + yn*l + zn*l2;

    
		for(j=0;j<t_neigh;++j)
		{
			if (s[i]!=s[neigh[j]])
			{
				num_interface++;

				switch (s[i]+s[neigh[j]])
				{
				// s=1 (agua) -- s=0 (ar) -- s=2 (oleo)  --- s=9 (solido)
				case 1: //agua - ar
				total+=0.5*(eps_WG); 
				ene_interface+=0.5*(eps_WG);
				break;
				case 2: //oleo -ar
				total+=0.5*(eps_OG);
				ene_interface+=0.5*(eps_OG);
				break;
				case 3: //agua - oleo
				total+=0.5*(eps_WO);
				ene_interface+=0.5*(eps_WO);
				break;
				case 9: //ar - solido
				total+=0.5*(eps_SG);
				ene_interface+=0.5*(eps_SG);
				break;
				case 10: //agua - solido
				total+=0.5*(eps_SW);
				ene_interface+=0.5*(eps_SW);
				break;
				case 11: //solido - oleo
				total+=0.5*(eps_SO);
				ene_interface+=0.5*(eps_SO);
				break;
				} //fim do SWITCH

			}//fim do IF

		}//fim do FOR


		if (s[i]==1) 
		{
			total+=Gw*(i/l2);
			ene_grav+=Gw*(i/l2);
		}//fim do if

		if (s[i]==2) 
		{
			total+=Go*(i/l2);
			ene_grav+=Go*(i/l2);
		}//fim do if


	} //fim do FOR

	total+=Lambda_w*(vol-t_vol)*(vol-t_vol);
	total+=Lambda_o*(vol_o-t_vol_o)*(vol_o-t_vol_o);

	fprintf(fout,"Total          : %f \n",total);
 	fprintf(fout,"Interfacial    : %f        num_interface: %ld\n",ene_interface,num_interface/2/t_neigh );
	fprintf(fout,"Gravitacional  : %f             int/grav: %f\n",ene_grav,ene_interface/ene_grav);
	fprintf(fout,"Lagrange Total : %d \n",Lambda_w*(vol-t_vol)*(vol-t_vol));
	fprintf(fout,"Lagrange Fração: %d \n",Lambda_o*(vol_o-t_vol_o)*(vol_o-t_vol_o));
	fprintf(fout,"\n\n\n");

	return total;
}

/*********************************************************************************************
*                       Rotina que escreve a conf no formato original .dsf                   *
*********************************************************************************************/
void save_conf(int num_steps,int iout) 
{
  
	int i;
  
	if(iout==1) 
	{	
		sprintf(output_file1,"%sgota_3d_L_%d_R_%d_a_%d_h_%d_w_%d_fo_%3.2f_LAST.dsf",CI,l,rg,a,h,w,Aw);
		flast = fopen(output_file1,"w");
		fprintf(flast,"# %d\n",num_steps);


		for(i=0;i<l3;i++) 
		{
			if(s[i]==1 || s[i]==2 ) 
			{
				fprintf(flast,"%d  %d\n",i, s[i]);
			}//fim do IF
		} //fim do FOR

  		fflush(flast);
    
  		fclose(flast);
	} //fim do if

	else if(iout==0) 
	{	

		fprintf(fconf,"# tempo  %d\n",num_steps);
		for(i=0;i<int_label;i++) 
		{
			if(s[w_inter[i]]==1 || s[w_inter[i]]==2 ) 
			{
      			fprintf(fconf,"%d  %d\n",w_inter[i], s[w_inter[i]] );

  			} //fim o IF

		} //fim do FOR

  		fflush(fconf);
  		fprintf(fconf,"\n\n"); 

	} //fim do ELSE IF

	return;
}

/*********************************************************************************************
*                              Rotina que testa a energia do flip                            *
*********************************************************************************************/
void teste_flip_droplet(int site, int s_teste, double delta_s, double delta)
{

	int nw,no,ng,ns;
	double ei, ef;
	
	nw = nv_w[site];
	ng = nv_g[site];
	no = nv_o[site];
	ns = nv_s[site];


	if(s[site]==1) ei= ng*eps_WG + ns*eps_SW + no*eps_WO;
	if(s[site]==0) ei= nw*eps_WG + ns*eps_SG + no*eps_OG;

	ef = nw*eps_WO + ng*eps_OG + ns*eps_SO;

	fprintf(foil ,"%10d %2d %2d %3d %3d %3d %3d %6.2f %6.2f %6.2f %6.2f %6.2f\n", site, s[site], s_teste, nw, no, ng, ns, ei, ef, ef-ei,delta_s, delta);

	return;
}

/*********************************************************************************************
*                              Rotina que testa a energia do flip                            *
*********************************************************************************************/
void teste_flip_surface(int site, int s_teste, double delta_s, double delta)
{

	int nw,no,ng,ns;
	double ei, ef;
	
	nw = nv_w[site];
	ng = nv_g[site];
	no = nv_o[site];
	ns = nv_s[site];


	if(s_teste==1) ef= ng*eps_WG + ns*eps_SW + no*eps_WO;
	if(s_teste==0) ef= nw*eps_WG + ns*eps_SG + no*eps_OG;

	ei = nw*eps_WO + ng*eps_OG + ns*eps_SO;

	fprintf(fwat ,"%10d %2d %2d %3d %3d %3d %3d %6.2f %6.2f %6.2f %6.2f %6.2f\n", site, s[site], s_teste, nw, no, ng, ns, ei, ef, ef-ei, delta_s, delta);

	return;
}











