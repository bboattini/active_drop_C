#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)

c=['#000000','#0000CD','#00FFFF','#008B00','#FFD700','#FF7F00','#CD0000']

#=========================================================================================
#=                                   Lendo o arquivo									 =
#=========================================================================================

arquivo = 'file.dsf'

t, V, Vw, Vo, E, bxw, byw, bxo,byo, rxw, ryw, rxo, ryo, txw, tyw,txo, tyo, npxw,npyw, npxo,npyo, vbw, vbo, vrw,vro, pw, po, fw, fo, vpw, vpo = np.loadtxt(arquivo, unpack=True)

#=========================================================================================
#=                                   Gráfico Volumes									 =
#=========================================================================================

grafico = 'volume_total.eps'

plt.plot(t,V,'-k', lw=1, label='Total V')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'$V$',fontsize=26)

plt.legend()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                   Gráfico Energia									 =
#=========================================================================================

grafico = 'energia.eps'

plt.plot(t,E,'-k', lw=1)

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'$E$',fontsize=26)

plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                     Raio da Base			    						 =
#=========================================================================================

grafico = 'B.eps'

plt.plot(t,bxw,'-' , c='tab:blue', lw=1, label='Bx Water')
plt.plot(t,byw,'--', c='tab:blue', lw=1,label='By Water')
plt.plot(t,bxo,'-' , c='tab:orange', lw=1, label='Bx Oil')
plt.plot(t,byo,'--', c='tab:orange', lw=1,label='By Oil')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'B',fontsize=26)

plt.legend()
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                   Gráfico raios gota									 =
#=========================================================================================

grafico = 'raio.eps'

plt.plot(t,rxw,'-' , c='tab:blue', lw=1, label='Rx Water')
plt.plot(t,ryw,'--', c='tab:blue', lw=1,label='Ry Water')
plt.plot(t,rxo,'-' , c='tab:orange', lw=1, label='Rx Oil')
plt.plot(t,ryo,'--', c='tab:orange', lw=1,label='Ry Oil')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'$R$',fontsize=26)

plt.legend()
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                   Gráfico Angulos									 =
#=========================================================================================

grafico = 'theta.eps'

plt.plot(t,txw,'-' , c='tab:blue', lw=1, label=r'$\theta _x$ Water')
plt.plot(t,tyw,'--', c='tab:blue', lw=1,label=r'$\theta _y$ Water')
plt.plot(t,txo,'-' , c='tab:orange', lw=1, label=r'$\theta _x$ Oil')
plt.plot(t,tyo,'--', c='tab:orange', lw=1,label=r'$\theta _y$ Oil')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'$\theta$',fontsize=26)

plt.legend()
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                    Gráfico vb    									 =
#=========================================================================================

grafico = 'v_baixo.eps'

plt.plot(t,vbw,'-' , c='tab:blue', lw=1, label=r'$V _{pil}$ Water')
plt.plot(t,vbo,'-' , c='tab:orange', lw=1, label=r'$V _{pil}$ Oil')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'V$_{p}$',fontsize=26)

plt.legend()
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

#=========================================================================================
#=                                    Gráfico vb    									 =
#=========================================================================================

grafico = 'v_per.eps'

plt.plot(t,vpw,'-' , c='tab:blue', lw=1, label=r'$V _p$ Water')
plt.plot(t,vpo,'-' , c='tab:orange', lw=1, label=r'$V _p$ Oil')

plt.xlabel(r'$t$',fontsize=26)
plt.ylabel(r'V$_p$',fontsize=26)

plt.legend()
plt.savefig(grafico, format='eps', dpi=400, bbox_inches='tight')
	
plt.close()

