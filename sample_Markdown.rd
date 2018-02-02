---
title: Test RMarkdown
date: 5 Jan 2017
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    code_folding: hide
    theme: cosmo
---

```{r global_options, include = FALSE}
knitr::opts_chunk$set(	fig.width	= 10, 
						fig.height	= 7, 
						fig.path	= "/Users/aesin/Desktop", 
						fig.align	= 'center', 
						dpi			= 300, 
						cache.path	= "/Users/aesin/Desktop", 
						warning		= TRUE, 
						message		= TRUE,
						tidy		= TRUE)
```

# Main Header 1 {.tabset .tabset-fade}
## Second header - tabbed
Here is some text that will appear in HTML. I can make this _bold_ or __italic__
Below is a code chunk. Cache is FALSE here, but typically should be TRUE

```{r nameOfMyChunk, warning = FALSE, message = FALSE, cache = FALSE}
x <- rnorm(100)
plot(density(x))
```

```{r functions, warning = FALSE, message = FALSE}
library(ape)
```

## Third header - tabbed