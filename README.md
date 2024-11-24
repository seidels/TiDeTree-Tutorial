---
author: Sophie Seidel
level: Intermediate
title: TiDeTree Tutorial
subtitle: Reconstruction time-scaled single-cell phylogenies from genetic lineage tracing data
beastversion: >= 2.7.
---


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

- BEAST 2 version 2.7 or later. [http://www.beast2.org/](http://www.beast2.org/)
- Tracer version 1.7 or later. [https://github.com/beast-dev/tracer/releases/latest](https://github.com/beast-dev/tracer/releases/latest)

# TiDeTree Package Installation [TODO]

# Setting up the analysis


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

