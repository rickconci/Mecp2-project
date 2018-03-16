# Mecp2-project
Part II project investigating the development of neural networks in dissociated cortical cells from Mecp2-deficient and wild-type mice. Recording spontaneous neuronal activity with multi-electrode arrays to compare functional connectivity and network topologies utilising Graph Theory. Coding in Matlab and R.  Supervised by Dr. Susanna Mierau &amp; Professor Ole Paulsen and collaborated with Manuel Schroeter from Bullmore Lab.

# Background 
Rett syndrome is an autism spectrum disorder caused by mutation in the MECP2 gene
![Imgur](https://i.imgur.com/7ga5zNh.png)

Understanding how network development is disrupted in Rett syndrome is key to finding and testing new treatments
![Imgur](https://i.imgur.com/ngr1on4.png)

# Aims
1. Investigate how networks develop and how this is disrupted in Rett Syndrome
2. Apply graph theory to investigate network properties 
3. Create model to distinguish  wild type and Mecp2-deficient networks

# Methods
We record 12 minutes worth of spontaneous activity from each microelectrode array once a week 
![Imgur](https://i.imgur.com/b49yN7r.png)

# Spike Detection 
Work done by @Timothysit (find code on Github)

Spike detection methods involve applying a band pass filter to the raw data and setting a threshold. We are currently developing a novel approach for spike detection based on Gaussian Mixture Models.
![Imgur](https://i.imgur.com/2VUUVvm.png)

By overlaying detected spikes on a picture of the neurons, we can eliminate noisy electrodes where there are no cells. 
![Imgur](https://i.imgur.com/uzh15Ea.png)

# Network Analysis
![Imgur](https://i.imgur.com/wMWD9DN.png)

The spiking activity can be visualised with a Raster Plot, with electrode number on the Y axis and time on the X axis. Each dark line denotes a certain number of spikes per time bin.  From this plot it is possible to visually identify possible network bursts (red rectangle).  Identifying which electrodes pick up bursting activity gives a first hint at the functional connections between the neurons on the MEA. 

### For Burst detection I used the ["Parameters for Burst Detection" by Bakkum et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3915237/pdf/fncom-07-00193.pdf) 
- Turn the binary spike activity into a vector of spikes times
- Create a histogram of the inter-spike-intervals for a specific number of spikes N
- Depending on the distribution of the ISI on the histogram select a threshold 
- Using the threshold we locate the start of each burst and save it as a separate spike train

### Correlated activity during a network burst as a sign of functional connectivity 

To consider the relationship between each region of the MEA we obtain the Pearson Correlation Coefficient between each pair of electrode recordings. This correlation is visualised by an Association Matrix. By considering an electrode as a node and a correlation with another electrode as an edge, we can create a weighted functional connectivity network. To make sure the edges  are signifying more trustworthy biological connections, we only select correlated nodes that are statistically significant. 

In the diagram below by [Schroeter et al 2015](http://www.jneurosci.org/content/35/14/5459.full.pdf), the same methods were used to investigate the development of mice hippocampal networks over time. 

![Imgur](https://i.imgur.com/mR132r7.png)

### Different measures of centrality allow us to characterize network hubs by function

![Imgur](https://i.imgur.com/YvJsGFP.png)
![Imgur](https://i.imgur.com/kFJ55Ox.png)


### Degree distributions define network topology and hubs
The degree, ki , of node i is the number of unique edges connecting node I with all other j = 1 …N-1 nodes. The degree distribution shows us the probability of a node in the network having a certain degree. In random networks, this probability distribution can be modelled by a Poisson distribution. In scale-free networks, such as the world wide web, the distribution follows a power law. Most biological networks have a broad scale distribution: at low degrees the distribution follow a scale-free  distribution and at high degrees  it switches to a small world distribution. Figure adapted from  Vella et al 2017. 

![Imgur](https://i.imgur.com/dFw6shE.gif)


# Future Directions

### Single cell resolution
Shifting from a 60 MEA to a 4096 electrode MEA would give us single cell resolution in our recordings and a much better signal to noise ratio.  With improved resolution we would be able to create biologically relevant weighted and directed functional networks.

![Imgur](https://i.imgur.com/2qVY6y9.png)


### Cell-type resolution
We could use optogenetics to selectively activate or inhibit the PV-ArchT inhibitory interneurons in our network. This approach allows us to alter the excitatory-inhibitory balance of the network.

