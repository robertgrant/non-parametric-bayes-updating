# non-parametric-bayes-updating
Non-parametric kernel methods for updating a Bayesian model with frequent batches or windows of data

©BayesCamp Ltd, CC-BY-NC-SA

This repository contains digital material associated with a paper currently submitted. The journal of choice does not permit pre-print archiving, but a post-print will appear here if it is accepted. You can also write to me for a copy.

## Abstract

One of the big attractions for people adopting Bayesian methods is the promise of "updating" their parameter estimates and predictions as more data arrive. Yesterday's posterior becomes today's prior. In practice, this is not simple, requiring at the very least a complete set of sufficient statistics, random samples from an unchanging population, and no changes of probability distribution for the prior.
In this day and age, there is great interest in big data (including fast, streaming data) but complex Bayesian models are usually prohibitively slow as the data accumulate. Another trend is to use windowing to privilege more recent data; there is no practicable Bayesian methodology for this as part of updating.
Sometimes, one would like to update without imposing an a priori distribution on yesterday's posterior and without estimating lots of statistics. This also might allow for faster analyses as new data arrive.
I discuss a kernel density approach in MCMC and SIR, consider two ways to reduce computational burden, evaluate different bandwidth strategies, compare them in two simulation studies and explore their feasibility in a real-life application.
The approach shows some promise. SIR is very fast but prone to problems of degeneracy and bandwidth shrinking. However, with human monitoring and adjustment, these can potentially be overcome. I identify some priorities for further research in the paper.


## Kernels and sampling algorithms

There is an existing and obvious parametric updating algorithm. Data batch 1 produces posterior sample 1, and from that we obtain summary statistics that get plugged into the parameters of the prior for data batch 2 (possibly via method of moments).

We want a different, non-parametric, approach, that is:
* flexible with regard to distributions
* fast
* and accurate.
We might not get all three, but we want to understand where the compromises arise and how to choose which aspect to prioritise. Our method will involve taking the posterior sample from batch b-1 (each "sample" comprises a number of "draws" from an unknown multivariate distribution) and using it as a prior sample for batch b.

This project tests out some kernel density approaches.

First, evaluating prior density by adding and normalising unbounded kernel densities (e.g. multivariate normal: MVN) around each of the prior draws, and adding a kernel density estimate, evaluated over the draws of yesterday's posterior, of the log-prior to today's log-likelihood (Algorithm 2).

But maybe we can make this run more quickly. Maybe we don't have to evaluate every prior draw's kernel. Two ideas came to mind:
* dividing the draws of yesterday's posterior into hypercuboids of parameter space, and calculating the kernel from the centroids of those, weighted by the number of draws in each (Algorithm 3)
* using a bounded kernel (e.g. a translated and rescaled beta distribution), so that only "neighbour" prior draws close to the current proposal (within each Metropolis-Hastings iteration) need to be evaluated (Algorithm 4)

Unfortunately, these have serious limitations. The number of hypercuboids will soon exceed the number of prior draws, and unless there is a fine grid, it introduces a probably unacceptable amount of distortion. To use bounded kernels, we have to find neighbours, which is computationally expensive in itself.

Apart from MCMC, we can also use sampling importance resampling (SIR). Take the prior draws, evaluate the likelihood at each draw, and store a number of posterior draws, proportional to the likelihood, at or around each one. We can apply the kernel for smoothing either before (SIR-B, smoothing the prior) or after (SIR-A, smoothing the posterior) (Algorithm 5)

Unfortunately, high-dimensional parameter spaces create two challenges. Firstly, as dimension increases, it is more likely that only one prior draw will produce posterior draws (degeneracy). This is a problem because you just get the kernel back, centred on that one draw, introducing bias in the location and spread, and loss of correlation. A sensible way around this is to "locally oversample" in the region of the MLE. An easy way to do this is to draw from truncated kernels inside a hypercuboid defined by, for example, the MLE +/- 3 standard errors (these being estimated from a frequentist analysis of batch b alone).

However, the second problem is that, as dimension increases, a larger proportion of the hypercuboid will lie outside a convex likelihood surface, and most plausible non-convex ones too. So, we may need some method to establish a more flexible region via likelihood hulls. This is not the focus of this project, but a sketch algorithm can be worked out (SIR-hulls-B, Algorithm 6, not included in the paper).


## Windowing

Suppose we have updated over batches 1 to 10. Now, we want to update our posteriors to include likelihood contributions from batches 11 onwards, but only to use the ten most recent batches each time. There are two ways to do this.

The first is simple, but needs more hardware and storage. It is to keep a collection of sliding windows (see https://www.oreilly.com/radar/the-world-beyond-batch-streaming-101/). In this example, you would keep updates from batch i to batch i+9. At any given batch b, you need to update i where b-9<=i<=b.

The second is more theoretically pleasing but at risk of distortion in real life. When you get to batch 11, you add its log-likelihood, then subtract batch 1's log-posterior, then add the original log-prior (batch 0, if you like) back in. You can even imagine having some kind of exponentially weighted likelihood.

## Bandwidth

I used Sheather-Jones bandwidths in each dimension in this project. I also tried out multiplying them by various constants. Sometimes, multimodality sneaks through, so I wrote some code to look for it using Hartigans' dip test and to inflate the bandwidth if that came up p<0.05. Some people might get agitated that this is not Bayesian, and likewise that simulation studies' long-run properties is somehow incompatible with Bayes. I direct you to http://www.robertgrantstats.co.uk/blog/26.html for my philosophical standpoint on probability, if you fit this profile.

If you want to make a contribution to this idea, probably the lowest-hanging fruit is to try out a local bandwidth algorithm like that proposed by Leo Breiman and co-authors way back in 1977.

## Method

I investigated these bandwidths, kernels, and samplers in three ways:

### prior-only simulation

I generated draws from three bivariate distributions (MVN, banana, multimodal typical of LASSO and other regularized models), and fed these in as prior sample 1. Then I updated ten times without any likelihood and compared distributions qualitatively (visually) and by 20-sliced 2-Wasserstein distances.
You can find this in the 1-prior-only folder. There is a subfolder for each of the bandwidths, and inside each, prior-only.R does the work and wassersteins.R compares the results.


### simulation study (logistic regression)

I ran a simulation study with a simple logistic regression scenario: 400 simulations, keeping 4000 draws at each update (and other keeping 2000 draws), using Algorithms 1-5.
You can find this in the 2-simulation folder. The simulation data are in data_batches and the first posterior for each simulation (from Stan) in first_posteriors. non-para-kernel-updating-sims.R does the sims (I actually set up four instances of R on a DigitalOcean Linux VM and ran this four times in parallel), recombine-simulations.R does what it says, along with generate-outputs.R and tables-figures.R. updates* image files show craft charts of the kernel densities of posteriors from each update of the first 30 simulations in different settings; the dark density superimposed on the tenth update is the Stan analysis of all data together, which we can regard as ground truth for comparison. means* image files show all 400 simulation means for beta0 changing with the updates. stats* image files show the distribution of posterior means across the simulations (mean +/- 2 SDs) for the three parameters over updates.


### feasibility case study (New York taxis)
As a feasibility case study, I applied SIR-A to Chris Whong's 2013 New York yellow taxi data, fitting a 50-parameter model to every taxi journey each day from 1--14 January 2013. When we hit degeneracy, I tried local oversampling (SIR-B). 1--7 Jan were updated.
You can find this in the 3-NYCtaxis folder. taxi_data* CSV files are the data in handy cleaned and reduced format. NYC-taxi-analysis.R does the work an is the place to go for readily adaptable code for SIR-A/B. data-prep.R produces the taxi_data* files from Chris Whong's files, so it is here just for completeness and audit trail. summary* image files show the posterior mean +/- 2 posterior SDs for each parameter in turn ober the updates. mean_evolution and SD_evolution show all the parameters' posterior stats superimposed over the updates. bivariate* takes a pair of parameters from the areas, the weekdays and the hours, and shows the posterior draws after each update. The red polygons show the convex hull of the following update, so you can see how it contracts over time. Problems evident in this are discussed in the paper.


## Conclusions

Kernel updating is a potentially useful trick when we need to scale up a serious Bayesian analysis to large datasets and run it regularly, as is often the case in commercial data science settings, not to mention pandemics of course.

It does not seem to have been done before; the only similar work I found was by Neufeld and Deutsch at the University of Alberta (see the 2006 paper "Data Integration with Non-Parametric Bayesian Updating"), although they used kriging on quantiles of univariate CDFs rather than kernels on posterior draws.

The submitted paper sets out more detailed recommendations and ares for future research.

You are welcome to get in touch (BayesCamp Ltd, 80 High Street (Unit 83), Winchester, SO23 9AT, England) to discuss this if it would be useful in your own work. But note:
* I do not have students
* I do not join bids for academic research funding
* I do not work for free
* I do not want a job
