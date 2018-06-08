# Quantum Monte Carlo for Interacting Fermions

This is a MSc Physics Engineering final project consisting of applying Quantum Monte Carlo (QMC) to a problem of interacting fermions.

Capturing the effects of electron correlations is not an easy task. The difficulty lies in devising a numerical method to solve the many-body Schrödinger equation in a reasonable amount of computer time. Naive methods have exponential complexity in the system size. The development and application of unbiased methods is a central point in correlated electron systems, particularly in 2D.

QMC methods are among the few unbiased methods available to date. In general, they circumvent the exponential complexity hurdle making it algebraic instead. However, for fermionic systems, a sign oscillation deems the algorithm exponential again and hence not very effective. This is due to the antisymmetric nature of the many-fermion wavefunction. One of the main tasks of this project is to investigate how to solve this issue. Some approximations exist and they perform differently depending on the problem at hand.

This masters thesis is about implementing an algorithm to deal with tight-binding problems for 2D interacting electronic models. Ultimately, the goal is to write a code that simulates a specific interacting electron system: a transition metal dichalcogenide (TMD) nanoribbon.

TMD's are graphene-like 2D materials that are promising from both a theoretical and an application perspective. A nanoribbon is a 2D  nanostructure that is much longer on one direction than on the other (like a ribbon), so that electronic edge states become relevant and lead to unusual properties. In practice, this means that we can use periodic boundary conditions on the longer direction, and open boundary conditions on the other.

An example of an interesting property of these nanostructures that we wish to investigate is magnetism. Furthermore, one might be interested in the different phases that arise within the system and in how do the transitions between them occur. For example, recent papers point at the possibility of topological superconductivity in TMD nanoribbons. This is a many-body effect that is only captured numerically by state of the art techniques such as QMC.

In short, our aim is to carry out a theoretical study (with particular emphasis on numerical aspects) of the properties of a TMD nanoribbon - a graphene-like 2D nanostructure - where electron interactions are particularly relevant using a QMC method.

## Getting started

In an attempt to mimic an usual practice in computer science, I decided to write a readme file to help me (and others) navigate through my project folder. The idea is to organize and document the material I found and created over the course of my research on this topic. The content of the project folder is exhaustively explained and every relevant folder or file is briefly commented below.

### Thesis

This folder contains a template for the thesis struture as specified by Instituto Superior Técnico (IST) Lisboa. I modified the template and added my own content. It was modularized from the start so as to allow modifications on particular sections of the document individually. It complies with the latest available specifications by IST (2015/2016 academic year).

### Introductory reading material (by folders)

#### _reading_

##### qmc


**M. E. J. Newman and G. T. Barkema** “*Monte Carlo Methods in Statistical Physics*” : general discussion on MC methods.

**Tao Pang** “*An Introduction to Quantum Monte Carlo Methods*” : overview of QMC methods.

**Zhaojun Bai, Wenbin Chen, Richard Scalettar, Ichitaro Yamazaki** “*Numerical Methods for Quantum Monte Carlo Simulations of the Hubbard Model*” : this is a core text that explores the usage of QMC to treat the Hubbard model.

##### _strongly_correlated_electrons_numerics_


***David Sénéchal, André-Marie, Tremblay Claude Bourbonnais (editors)* , Zhang** (who wrote the chapter on QMC) “*Theoretical Methods for Strongly Correlated Electrons*” : proceedings of a 1999 workshop on Theoretical methods for strongly correlated electrons; a very general discussion of different methods, which might be useful to acquire a general perspective. Includes a chapter on auxiliary-field fixed node QMC (which overcomes the fermion sign problem).

## Built with

*LaTeX*

*C++*

*Python*

## Contributing

*Francisco Monteiro de Oliveira Brito* (CeFEMA, Instituto Superior Técnico de Lisboa, Centro de Física da Universidade do Porto)

## Versioning

v1 - latest update 23.02.2018

## Authors

*Francisco Monteiro de Oliveira Brito* (CeFEMA, Instituto Superior Técnico de Lisboa, Centro de Física da Universidade do Porto)

## Acknowledgements

*Miguel Amador* - thesis template
