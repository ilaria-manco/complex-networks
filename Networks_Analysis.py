"""
Networks_Analysis.py
Module to perform analysis on the data for the Networks Project.
I Manco    26/03/17

Classes: - Ensemble
         - Analysis

Comments:
The module log_bin_CN_2016.py was written by James Clough (jrc309@ic.ac.uk)
for the 2016 Complexity & Networks course

"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import log_bin_CN_2016 as lb
import BAmodel
reload(BAmodel)
        
class Ensemble:
    def __init__(self, a, k, max_k, N, m, attachment, runs, L = 0):
        """Initialise Ensemble.
        N          -- number of final nodes
        m          -- number of edges to add to new vertice
        attachment -- type of attachment 
        runs       -- number of graphs in the ensemble
        L          -- length of random walk
        """
        self.N = N
        self.m = m
        self.attachment = attachment
        self.L = L
        self.k = k
        self.max_k = max_k
        self.runs = runs
        self.a = a
        self.pdf = self.degree_prob()[0]
        self.cdf = self.degree_prob()[1]
        self.log_binned = self.log_bin(self.a)
        
    def degree_prob(self): 
        k = self.k
        pdf = np.bincount(k)/(float(self.N)*self.runs)
        cdf = []
        for degree in range(len(pdf)):
            cdf.append(sum(pdf[degree:]))
        return pdf, cdf   
    
    def log_bin(self, a):
        #x = self.k[5:]
        x = self.k
        vals, counts = lb.lin_bin(x, int(max(x)))
        b, c = lb.log_bin(x, 1., 1.5, a, debug_mode=False)
        slope = stats.linregress(np.log(b[4:]), np.log(c[4:]))
        #print np.log(b), np.log(c)
        return slope, b[1:], c[1:]
        #return slope, b, c
        
class Analysis:
    def __init__(self, attachment, analysis, L = 0):
        #self.m = m
        self.attachment = attachment
        self.analysis = analysis
        self.L = L
        self.ensembles = []
        if analysis == 'm':
            for n in range(6):
                runs = 2**(10-n)
                if self.attachment == 'preferential':
                    self.ensembles.append(BAmodel.get_ensemble(2**n, 1e5, 'preferential', self.analysis, runs))
                if self.attachment == 'random':
                    self.ensembles.append(BAmodel.get_ensemble(2**n, 1e4, 'random', self.analysis, runs))
        if analysis == 'finite size':
            #for n in range(2, 6):
            for n in range(2, 4):
                runs = 2**(13-n)
                #runs = 2**(8-n)
                N = 10**n
                if self.attachment == 'preferential':
                    self.ensembles.append(BAmodel.get_ensemble(4, N, 'preferential', 'finite size', runs))
                if self.attachment == 'random':
                    self.ensembles.append(BAmodel.get_ensemble(4, N, 'random', 'finite size', runs))
        if analysis == 'none':
            self.ensembles.append(BAmodel.get_ensemble(3, 1e4, 'random walk', 'none', 100, L = 0))
            self.ensembles.append(BAmodel.get_ensemble(2, 1e4, 'random walk', 'none', 100, L = 1))
            self.ensembles.append(BAmodel.get_ensemble(2, 1e4, 'random walk', 'none', 10, L = 10))
        
                
    def plot_max_k(self):
        N = []
        k_1 = []
        k_theor = []
        for ensemble in self.ensembles:
            N.append(ensemble.N)
            k_1.append(np.mean(ensemble.max_k))
            if self.attachment == 'preferential':
                theoretical = (1./2.)*(-1 + np.sqrt(1+4*ensemble.N*ensemble.m*(ensemble.m+1)))
            if self.attachment == 'random':
                m = ensemble.m
                theoretical = ensemble.m - (np.log(ensemble.N)/(np.log(m)-np.log(m+1)))
            k_theor.append(theoretical)
        plt.loglog(N, k_theor, 'b')
        fit = stats.linregress(np.log(N), np.log(k_1))
        plt.loglog(N, k_1, '--.')
        plt.xlabel('$N$')
        plt.ylabel('$k_1$')
        plt.show()
        return fit
    
    def plot_pdf(self):
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'black']
        col = 0
        for ensemble in self.ensembles:
            pdf_theor = []
            if ensemble.attachment == 'random walk':
                plt.loglog(ensemble.pdf,'.', label = 'L = '+ str(ensemble.L), color = colors[col])
            else:
                plt.loglog(ensemble.pdf,'.', label = 'm = '+ str(ensemble.m), color = colors[col])
            for k in range(int(np.min(ensemble.k)), int(1e4)):
                if ensemble.attachment == 'preferential':
                    pdf_theor.append((2.*ensemble.m*(ensemble.m+1))/float(k*(k+1)*(k+2)))
                if ensemble.attachment == 'random':
                    pdf_theor.append((1./(1.+ensemble.m))*((ensemble.m/(ensemble.m+1.))**(k-ensemble.m)))
            plt.loglog(range(int(np.min(ensemble.k)), int(1e4)), pdf_theor, '--', linewidth = 0.5, color = colors[col])
            col = col + 1
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.legend(prop={'size':14})
        plt.show()  
        
    def plot_binned_pdf(self):
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'black']
        col = 0
        for ensemble in self.ensembles:
            log_binned = ensemble.log_binned[2]
            log_binned_k = ensemble.log_binned[1]
            pdf_theor = []
            for k in range(int(np.min(ensemble.k)), int(1e4)):
                if ensemble.attachment == 'preferential':
                    pdf_theor.append((2.*ensemble.m*(ensemble.m+1))/float(k*(k+1)*(k+2)))
                if ensemble.attachment == 'random':
                    pdf_theor.append((1./(1.+ensemble.m))*((ensemble.m/(ensemble.m+1.))**(k-ensemble.m)))
            if self.analysis == 'm':
                labels = 'm = '+ str(ensemble.m)
            if self.analysis == 'finite size':
                if ensemble.N == 1e2:
                    labels = 'N = $10^2$'
                if ensemble.N == 1e3:
                    labels = 'N = $10^3$'
                if ensemble.N == 1e4:
                    labels = 'N = $10^4$'
                if ensemble.N == 1e5:
                    labels = 'N = $10^5$'   
                if ensemble.N == 1e6:
                    labels = 'N = $10^6$'  
            plt.loglog(log_binned_k, log_binned, label = labels, color = colors[col])
            plt.loglog(range(int(np.min(ensemble.k)), int(1e4)), pdf_theor, '--', linewidth = 0.5, color = colors[col])
            col = col + 1
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.legend(prop={'size':14})
        plt.show()  

    def collapse(self):
        for ensemble in self.ensembles:
            log_binned = np.array(ensemble.log_binned[2])
            log_binned_k = np.array(ensemble.log_binned[1])
            if self.attachment == 'preferential':
                k1 = (1./2.)*(-1 + np.sqrt(1+4*ensemble.N*ensemble.m*(ensemble.m+1)))
            if self.attachment == 'random':
                m = ensemble.m
                k1 = ensemble.m - (np.log(ensemble.N)/(np.log(m)-np.log(m+1)))
            x = log_binned_k/k1
            p_inf = []
            for k in log_binned_k:
                if self.attachment == 'preferential':
                    p_inf.append((2.*ensemble.m*(ensemble.m +1))/(k*(k+1)*(k+2)))
                if self.attachment == 'random':
                    p_inf.append((1./(1.+ensemble.m))*((ensemble.m/(ensemble.m+1.))**(k-ensemble.m)))
            p_inf = np.array(p_inf)
            y = log_binned/p_inf
            if ensemble.N == 1e2:
                labels = 'N = $10^2$'
            if ensemble.N == 1e3:
                labels = 'N = $10^3$'
            if ensemble.N == 1e4:
                labels = 'N = $10^4$'
            if ensemble.N == 1e5:
                labels = 'N = $10^5$'   
            if ensemble.N == 1e6:
                labels = 'N = $10^6$'  
            plt.loglog(x, y, '-o', label = labels)
        plt.ylabel('$P(k)$ / $P_\infty$$(k)$')
        plt.xlabel('$k$ / $k_1$')
        plt.legend(prop={'size':14})
        plt.show()
        
    def plot_cdf(self):
        colors = ['b', 'g', 'r', 'c', 'm', 'y']
        col = 0
        for ensemble in self.ensembles:
            cdf_theor = []
            if self.analysis == 'finite size':
                if ensemble.N == 1e2:
                    labels = 'N = $10^2$'
                if ensemble.N == 1e3:
                    labels = 'N = $10^3$'
                if ensemble.N == 1e4:
                    labels = 'N = $10^4$'
                if ensemble.N == 1e5:
                    labels = 'N = $10^5$'   
                if ensemble.N == 1e6:
                    labels = 'N = $10^6$' 
            if self.analysis == 'm':
                labels = 'm = '+str(ensemble.m)
            plt.loglog(ensemble.cdf, '.', label = labels, color = colors[col])
            for k in range(int(np.min(ensemble.k)), int(1e4)):
                cdf_theor.append((1./(ensemble.m*np.log(ensemble.m/(ensemble.m+1))))*((ensemble.m/(ensemble.m+1.))**(k-ensemble.m+1)))
            plt.loglog(range(int(np.min(ensemble.k)), int(1e4)), cdf_theor, '--', linewidth = 0.5, color = colors[col])
            col = col + 1
        plt.xlabel('$k$')
        plt.ylabel('$P(k < k^*)$')
        plt.legend(prop={'size':14})
        plt.show()
        
    def ks(self):
        for ensemble in self.ensembles:
            pdf_data = ensemble.degree_prob()[0][ensemble.m:]
            pdf_theor = []
            m = float(ensemble.m)
            for degree in range(1, max(ensemble.k)+1):
                if ensemble.attachment == 'preferential':
                    pdf_theor.append((2*m*(m+1.))/float(degree**3))
                if ensemble.attachment == 'random':
                    pdf_theor.append((1./(1.+ensemble.m))*((ensemble.m/(ensemble.m+1.))**(degree-ensemble.m)))
            ks = stats.ks_2samp(pdf_data[:10], pdf_theor[:10])
            pdf_theor = pdf_theor[ensemble.m:]
            p_value = []
            D = []
            degree = []
            for k in range(1, np.max(ensemble.k)+1):
                degree.append(k)
                p_value.append(stats.ks_2samp(pdf_data[:k+1], pdf_theor[:k+1])[1])
                D.append(stats.ks_2samp(pdf_data[:k+1], pdf_theor[:k+1])[0])
            plt.xlabel('$k$')
            plt.ylabel('$p-value$')
            plt.plot(degree, p_value, label = 'm = '+str(ensemble.m))
            plt.legend(prop={'size':14})
            print stats.ks_2samp(pdf_data[:10], pdf_theor[:10]), ks
