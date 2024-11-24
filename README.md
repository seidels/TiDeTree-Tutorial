---
author: Sophie Seidel
level: Intermediate
title: TiDeTree Tutorial
subtitle: Reconstruction time-scaled single-cell phylogenies from genetic lineage tracing data
beastversion: >= 2.7.
---


# Background

# Background

## Introduction

### Understanding cell dynamics

This tutorial focuses on using TiDeTree to address key questions in developmental biology and cellular population dynamics. Specifically, TiDeTree is designed to tackle the following challenges:

1. **Reconstructing Time-Scaled Single-Cell Phylogenies:** How can we accurately infer the evolutionary relationships among individual cells in a population using genetic lineage tracing data?
2. **Estimating Population Dynamics:** How can we quantify cell division, death, and differentiation rates from lineage tracing experiments to better understand cellular behaviors over time?

These questions are central to studying developmental processes, such as tissue formation, cancer progression, or immune cell dynamics.

## TiDeTree Foundations

TiDeTree is a BEAST 2 package designed for inferring time-scaled single-cell phylogenies and estimating population dynamics parameters such as cell division, death, and differentiation rates. It is specifically tailored for analyzing genetic lineage tracing data, where random edits introduced through technologies like CRISPR-Cas9 allow tracking of cell lineages over time. TiDeTree incorporates a specialized editing and silencing model to account for the unique features of lineage tracing data:

1. **Editing Events:** Genetic modifications (e.g., insertions or deletions introduced by CRISPR-Cas9) serve as markers of lineage history.

2. **Silencing Events:** The model accounts for the possibility that some edits may become undetectable due to biological processes like gene silencing.

This dual approach enables accurate reconstruction of cell lineages and robust estimation of underlying population dynamics.

## Objectives

This tutorial aims to equip users with the knowledge and skills to use TiDeTree for their own research. Specifically, by following the tutorial, users will:

1. Understand the purpose and applications of TiDeTree in single-cell phylogenetics and population dynamics.
2. Gain theoretical insights into Bayesian methods and the editing-silencing model used by TiDeTree.
3. Learn how to prepare genetic lineage tracing data for analysis.
4. Successfully install TiDeTree and set up the required software environment.
5. Run example analyses to infer phylogenies and estimate population parameters.
6. Interpret results and understand how to adjust parameter settings for different datasets.

This tutorial is designed for researchers working with single-cell data, developmental biology, or lineage tracing studies, and assumes basic familiarity with phylogenetic concepts and BEAST 2.

## Prerequisites

Before starting this tutorial, users should ensure they have:

- Basic knowledge of phylogenetics and Bayesian inference.
- Java 17 installed (required for BEAST 2).
- BEAST 2 and TiDeTree installed on their computer.
- Example genetic lineage tracing data formatted according to TiDeTree specifications.

# Programs used in this Exercise 

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014 --file Tutorial-Template/master-refs.bib %}. This tutorial uses the BEAST2 version 2.4.2.

### BEAUti – Bayesian Evolutionary Analysis Utility
BEAUti is a utility program with a graphical user interface for creating BEAST2 input files, which are written in XML. The eXtensible Markup Language (XML) is a general-purpose markup language, which allows for the combination of text and additional information. The use of the XML makes analysis specification very flexible and readable by both the program and people. The XML file specifies all the components of the analysis, including sequences, node calibrations, models, priors, output file names.


### Tracer
Tracer is used for assessing and summarizing the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and assessment of convergence and it also calculates 95% credible intervals (which approximate the 95% highest posterior density intervals) and effective sample sizes (ESS) of parameters. Contrary to the other software in this section, Tracer is not distributed with BEAST2 and needs to be downloaded separately [here](https://beast.community/tracer).
----

# Practical: Exercise title

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Fusce porta augue id vulputate iaculis. Sed non posuere lectus. Integer at magna quis nulla tempus cursus. Integer ut nisl elit. Nam pellentesque pharetra orci eu facilisis. Nullam vitae leo tempus nunc consequat finibus ut ut nisi. Phasellus vitae faucibus dolor, id venenatis lacus. Sed eu lacus a nibh luctus semper. Proin non tellus odio. Duis elementum lorem eget nisl rhoncus feugiat. In hendrerit vehicula purus. Aliquam ornare libero quis tincidunt efficitur. Nam sapien augue, mattis nec hendrerit ut, commodo in nisi.

## This is a subsection
Quisque a dictum erat. Curabitur congue sapien sit amet pharetra pretium. Proin posuere euismod velit, eget faucibus ex varius id. Fusce sodales maximus malesuada. Mauris auctor dui in justo interdum egestas. Cras dapibus commodo nulla vitae congue. Vestibulum sit amet justo sit amet ex pretium bibendum. Donec ac mollis lorem, vel semper enim. Suspendisse sit amet auctor dui. Nullam ac efficitur mauris. Proin aliquam tincidunt felis nec semper. Vestibulum vestibulum, eros sit amet consectetur blandit, elit dolor posuere sem, a porta purus odio sit amet quam. Quisque dapibus erat sem, at vulputate libero dapibus sit amet. Mauris rhoncus odio nisl, nec interdum lacus consequat nec. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos.

Etiam tincidunt porttitor rutrum. Nulla facilisi. Mauris vehicula, justo ac ultricies tempus, quam erat hendrerit dui, vel pharetra sapien nibh vel ex. Sed molestie eu dui in laoreet. Pellentesque ultrices, orci vitae lacinia suscipit, erat sapien elementum ligula, sit amet viverra ante lorem eget elit. Cras euismod felis libero, pharetra lobortis arcu congue vehicula. Nullam posuere dapibus mauris, eget vulputate ligula auctor eget. Aenean et tempus est. Aliquam vehicula arcu vitae metus dictum viverra. Aliquam vitae purus mauris. Nullam interdum mauris eget sagittis consequat. Quisque in orci elementum, eleifend tortor eget, bibendum orci. Etiam aliquet dolor non neque semper fermentum. Praesent vitae venenatis mi, ut faucibus ligula. Phasellus vitae lorem neque. Interdum et malesuada fames ac ante ipsum primis in faucibus.

### This is a sub-subsection
Etiam posuere urna ut condimentum sagittis. Suspendisse posuere, ex nec eleifend fringilla, nisl augue posuere augue, elementum mollis justo felis sed purus. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Mauris efficitur eros ut turpis elementum vestibulum. Sed sit amet nisi at nunc luctus laoreet id ac enim. Aliquam elementum risus id urna dictum fringilla. Aenean lobortis, risus euismod molestie pulvinar, massa odio pharetra nulla, vitae facilisis neque magna sed lorem. Praesent ipsum enim, commodo ut pharetra in, sollicitudin ac massa. Donec et interdum mauris. Ut molestie, risus quis fermentum placerat, diam risus posuere nisi, eget viverra tortor neque ac sem. Donec viverra magna non dolor aliquam, in suscipit massa facilisis. Suspendisse congue arcu sed risus consectetur commodo. Aenean metus odio, volutpat at tincidunt id, ullamcorper in dui. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras ut sem in odio sodales iaculis non quis neque.

Quisque non mollis massa, nec eleifend dolor. Proin porta elit metus, a lobortis enim venenatis ac. Nam scelerisque consectetur mi et gravida. Vestibulum placerat, est vitae euismod finibus, purus nisl viverra quam, eget condimentum mauris magna vel nisl. Phasellus pretium vitae diam in volutpat. Cras gravida non quam ut consectetur. Vivamus congue vulputate lorem.

## This is another subsection
Praesent sodales est in tempor commodo. Suspendisse nulla metus, gravida eget malesuada vel, viverra eu felis. In vitae leo facilisis, ornare nunc nec, tempor tortor. Duis pretium mi eros, at consequat neque tincidunt eget. Mauris vestibulum venenatis arcu, eget lacinia arcu faucibus ut. Phasellus aliquam dui ipsum, a eleifend lacus fermentum at. Suspendisse congue orci quis ante consequat ornare. Integer a massa blandit, vestibulum eros ut, pulvinar augue. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos.

Quisque a urna a massa congue rhoncus. Donec bibendum tempus velit. Nam varius augue sit amet lacinia hendrerit. Proin tincidunt massa ut mi vestibulum placerat. Phasellus eget dui molestie, aliquet libero efficitur, vehicula ex. Pellentesque ultricies ante leo, eu lobortis odio convallis id. Donec vitae risus dui. Nulla orci velit, ultricies sed finibus quis, blandit quis arcu. Morbi non neque non odio rutrum condimentum. Vivamus libero metus, vehicula vitae elit ac, tincidunt pretium dui. Proin condimentum fringilla diam, blandit blandit nisl dapibus vel. Proin ante felis, accumsan eget ligula et, lobortis dictum nunc. Mauris a ante dignissim ipsum tincidunt tristique.

-------

# Tutorial style guide

## Text styling

This is how to write _italic text_.

This is how to write **bold text**.

This is how to write **_bold and italic text_**.

Do text superscripts like this 7^th, x^2y or  x^(2y + 3z).


## Lists

### Unnumbered lists

- Lorem ipsum dolor sit amet, consectetur adipiscing elit.
- Integer pharetra arcu ut nisl mollis ultricies.
	- Fusce nec tortor at enim cursus dictum.
	- Phasellus nec urna quis velit eleifend convallis sodales nec augue.
- In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
- Nam vitae turpis eu lacus imperdiet mollis id at augue.
- Sed sed turpis ac dolor mollis accumsan.


### Numbered lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	1. Fusce nec tortor at enim cursus dictum.
	2. Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.

### Mixed lists

1. Lorem ipsum dolor sit amet, consectetur adipiscing elit.
2. Integer pharetra arcu ut nisl mollis ultricies.
	* Fusce nec tortor at enim cursus dictum.
	* Phasellus nec urna quis velit eleifend convallis sodales nec augue.
1. In iaculis turpis in massa facilisis, quis ultricies nibh ultricies.
1. Nam vitae turpis eu lacus imperdiet mollis id at augue.
1. Sed sed turpis ac dolor mollis accumsan.


## Figures


<figure>
	<a id="fig:example1"></a>
	<img style="width:25%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 1: This figure is 25% of the page width.</figcaption>
</figure>


<figure>
	<a id="fig:example2"></a>
	<img style="width:10%;" src="figures/Logo_bw.png" alt="">
	<figcaption>Figure 2: This figure is only 10% of the page width.</figcaption>
</figure>



# Code

A bit of inline monospaced font can be made `like this`. Larger code blocks can be made by using the code environment:

Java:

```java
public class HelloWorld {

    public static void main(String[] args) {
        // Prints "Hello, World" to the terminal window.
        System.out.println("Hello, World");
    }

}
```

XML:

```xml
	<BirthDeathSkylineModel spec="BirthDeathSkylineModel" id="birthDeath" tree="@tree" contemp="true">
	      <parameter name="origin" id="origin" value ="100" lower="0."/>    
	      <parameter name="R0" id="R0" value="2" lower="0." dimension ="10"/>
	      <parameter name="becomeUninfectiousRate" id="becomeUninfectiousRate" value="1" lower="0." dimension ="10"/>
	      <parameter name="samplingProportion" id="samplingProportion" value="0."/>
	      <parameter name="rho" id="rho" value="1e-6" lower="0." upper="1."/>
	</BirthDeathSkylineModel>
```

R:

```R
	> myString <- "Hello, World!"
	> print (myString)
	[1] "Hello, World!"
```

# Equations

Inline equations: {% eqinline \dot{x} = \sigma(y-x) %}

Displayed equations: 
{% eq \left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right) %}



## Instruction boxes

Use block-quotes for step-by-step instruction that the user should perform (this will produce a framed box on the website):

> The data we have is not the data we want, and the data we need is not the data we have.
> 
> We can input **any** formatted text in here:
>
> - Even
> - Lists
>
> or equations:
>
> {% eq (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right) %}






# Hyperlinks

Add links to figures like this: 

- [Figure 1](#fig:example1) is 25% of the page width.
- [Figure 2](#fig:example2) is 10% of the page width. 

Add links to external URLs like [this](http://www.google.com). 

Links to equations or different sections within the same document are a little buggy.


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Tutorial-Template/master-refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- BEAST 1 website and documentation: [http://beast.bio.ed.ac.uk](http://beast.bio.ed.ac.uk)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users) 

----

# Relevant References

{% bibliography --cited --file Tutorial-Template/master-refs.bib %}

