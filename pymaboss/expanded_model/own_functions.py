# coding: utf-8

# Functions
import maboss
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def barplot_last_nodes_probtraj_selected(df,sel_nodes): 
    last_nodes=df.get_last_nodes_probtraj()
    truthvals=df.get_last_nodes_probtraj().columns.str.contains(str.join("|",sel_nodes))
    last_nodes=last_nodes.loc[:,truthvals]
    plt.barh(np.arange(len(last_nodes.columns)),last_nodes.values[0],tick_label=last_nodes.columns)
    return last_nodes

# any number of inputs
def fcn_truthval_states(df,string_cellfates):
    truth_vals=[0]*len(string_cellfates)
    for n in range(0,len(string_cellfates),1):
        truth_vals[n]=df.columns.str.contains(string_cellfates[n])
    return truth_vals

# plot cell fate dynamics
def plot_cell_fates(df,model_name,str_plotname,plot_pars,list_truthval,list_labels,save_flag):
    plt.figure(figsize=(plot_pars[0],plot_pars[1])); fontsize_val=plot_pars[2]; 
    plot_name=model_name+str_plotname; linewidth_val=plot_pars[3]
    for n in range(0,len(list_labels),1):
        plt.plot(np.sum(df.loc[:,list_truthval[n]],axis=1),label=list_labels[n],linewidth=linewidth_val)
    plt.xlabel("time",fontsize=fontsize_val); plt.ylabel("probability",fontsize=fontsize_val);
    plt.title(plot_name,fontsize=fontsize_val); 
    plt.legend(bbox_to_anchor=(1.04,0.5),loc="center left",borderaxespad=0,fontsize=fontsize_val)
    plt.ylim(0,plot_pars[4])
    if save_flag:
        plt.savefig("figures/"+plot_name+".png",bbox_inches='tight')
    if len(plot_pars)>5:
        plt.xlim(0,plot_pars[5])

def plot_node_dynamics(maboss_results_object,nodes_string,fig_pars): 
    nodes_probtraj_df=pd.DataFrame(maboss_results_object.get_nodes_probtraj())
    col_inds=nodes_probtraj_df.columns.str.contains('|'.join(nodes_string),case=False)
    nodes_probtraj_df.loc[:,col_inds].plot(linewidth=fig_pars[0])
    fig_nodetraj=plt.gcf(); fig_nodetraj.set_size_inches(fig_pars[1],fig_pars[2])
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.xlabel(fig_pars[3]); plt.ylabel(fig_pars[4]); 
    if len(fig_pars)>5:
        plt.xlim(0,fig_pars[5])

# plot states with prob larger than x
def plot_states_larger_prob_x(maboss_results_object,x,fig_pars):
    states_probtraj=maboss_results_object.get_states_probtraj()
    states_probtraj.loc[:,(np.sum(states_probtraj>x)>0) &                         (states_probtraj.columns!='<nil>')].plot(linewidth=fig_pars[0])
    fig_statetraj=plt.gcf(); fig_statetraj.set_size_inches(10,6);
    plt.legend(bbox_to_anchor=(1.04,0.5),loc="center left",borderaxespad=0)
    fig_statetraj=plt.gcf(); fig_statetraj.set_size_inches(fig_pars[1],fig_pars[2])
    plt.xlabel(fig_pars[3]); plt.ylabel(fig_pars[4])
    if len(fig_pars)>5:
        plt.xlim(0,fig_pars[5])

# solve ODE with init conds
def solve_ode(ode_object,Tmax,n_step,initvals):
    ode_object.set_initial_value(initvals, 0)
    x=[]; t=[]; dt=Tmax/n_step
    while ode_object.successful() and ode_object.t <= Tmax:
        ode_object.integrate(ode_object.t+dt)
        x.append(ode_object.y);  t.append(ode_object.t)
    return t,x

def fcn_sgn_eps_approx(x,epsilon):
    return x/np.sqrt(np.square(x) + epsilon**2)

def fcn_sgn_hvs_eps_approx(x,epsilon):
    return 0.5*fcn_sgn_eps_approx(x,epsilon)*(1 + fcn_sgn_eps_approx(x,epsilon))

def fcn_or_operator_fcns(f1,f2,epsilon):
    return fcn_sgn_eps_approx(f1+f2,epsilon)

# def f1_ode(t,u):
#    return [u[1], 3*(1 - u[0]*u[0])*u[1] - u[0]]

# defined ODE for cell fate
def f_cellfate_ode(x,t,params,ko_flag):
    kras=0; dnadam=1; dnarep=2; chek1=3; mitosis=4; celldeath=5;
    tau=params[0]; tau_decay=params[2];
    xdot=np.zeros(6);
    if params[1]==1: 
        xdot[kras]=-x[kras]/tau_decay
    else:
        xdot[kras]=(1-x[kras])/tau_decay
    epsilon=params[3]; threshold=params[4];
    # dnadam
    f_growth_dnadam=(1-x[dnarep])*x[kras];
    f_all_dnadam=(f_growth_dnadam - x[dnadam] )/tau;
    inh_terms=x[dnarep];
    s_h_f_all_dnadam=fcn_sgn_hvs_eps_approx(f_all_dnadam,epsilon);
    decrease_truth=fcn_or_operator_fcns(fcn_sgn_eps_approx(inh_terms,epsilon),s_h_f_all_dnadam,epsilon);
    increase_truth=fcn_or_operator_fcns(fcn_sgn_eps_approx(f_growth_dnadam,epsilon),1-s_h_f_all_dnadam,epsilon);
    xdot[dnadam]=f_all_dnadam*decrease_truth*increase_truth;
    # dnarep
    xdot[dnarep]=(x[dnadam]- x[dnarep])/tau;
    # chek1
    xdot[chek1]=(x[dnarep]- x[chek1])/tau;
    # mitosis
    growth_terms_mitosis = (1-x[chek1])*x[kras];
    f_all_mitosis=(growth_terms_mitosis - x[mitosis] )/tau;
    s_h_f_all_mitosis=fcn_sgn_hvs_eps_approx(f_all_mitosis,epsilon);
    decrease_truth_mitosis=s_h_f_all_mitosis;
    increase_truth_mitosis=fcn_or_operator_fcns(fcn_sgn_eps_approx(growth_terms_mitosis,epsilon),1-s_h_f_all_mitosis,epsilon);
    xdot[mitosis]=f_all_mitosis*decrease_truth_mitosis*increase_truth_mitosis;
    # celldeath=celldeath|(mitosis&dnadam)
    k_cd=params[5]; # k_cd for cell death rate=x/(x+k_cd)
    f_mitosis=(1+k_cd)*(x[mitosis])/(x[mitosis]+ k_cd);
    f_dnadam_celldeath=(1+k_cd)*(x[dnadam])/(x[dnadam] + k_cd);
    f_growth_celldeath=f_mitosis*f_dnadam_celldeath;
    f_all_celldeath=(f_growth_celldeath - x[celldeath])/tau;
    s_h_f_all_celldeath=fcn_sgn_hvs_eps_approx(f_all_celldeath,epsilon);
    increase_truth_celldeath=fcn_or_operator_fcns(fcn_sgn_eps_approx(f_growth_celldeath,epsilon),1-s_h_f_all_celldeath,epsilon);
    decrease_truth_celldeath=s_h_f_all_celldeath;
    xdot[celldeath]=f_all_celldeath*decrease_truth_celldeath*increase_truth_celldeath;
    if ko_flag:
        xdot[ko_flag]=0
    return xdot
    