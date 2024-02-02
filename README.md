# SpikeSorting

Development of an algorithm to do spike sorting and automatic curation of a signal recorded from the brainstem of mice. It is used to isolate single-unit activity in order to use it in medical technologies to decode signals from individual neurons.

# Motivation

Over the years, neuroscientists developed plenty of experimental methods to study how the brain works. One of the most used methods is extracellular recording. It can be used in
medical technologies as a treatment of disorders such as epilepsy and paralysis and for these treatments, single-unit activity is very important and of particular interest. This activity can be required by researchers in basic science to study how a type of neuron responds to a specific stimulus, while most medical technologies will need it for a “decoding” algorithm that operates on signals from individual neurons to decode movements or intentions.

The issue with this is that often, the recorded signal is a multi-unit activity, meaning that is it the sum of the signal from several (2 to 10) neurons surrounding the electrode. The reason behind this is the size of the recording probe and in such a case, the solution would be to use spike sorting which is the process of separating the recorded units into single units. Several algorithms for spike sorting have been published and are available online, but none is considered as a universally accepted solution.

This code can be used to differentiate, within one recording, different brain activities with high accuracy. 

