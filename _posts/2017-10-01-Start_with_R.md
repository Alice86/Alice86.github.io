R Set-up
--------

    # 安装 rticles 包 
    install.packages("devtools")
    library(devtools) 
    devtools::install_github("rstudio/rticles")
    install.packages("quantmod", repos = "http://R-Forge.R-project.org")

Some tools for objects managment
--------------------------------

    # 注意:R区分大小写，R 里面的目录要用反斜杠/或者\\ 
    getwd() 
    setwd("D:\\RPROJECT\\DataAnalysis\\data")
    ls()
    rm() 
    options(digits=3) 
    save.image("filename")

Source Codes
------------

View source codes is a powerful way of learning.

-   Simple functions: function.name, no ()  
-   Class Function: methods()

<!-- -->

    summary
    methods(summary)
    summary.lm
    # Non-visible functions are asterisked *
    getAnywhere()
    getS3method()

-   Download \*.tar.gz for source codes.  
    R main site: <http://www.r-project.org/>,  
    Download - CRAN - Source Code for all Platforms - sorted by name -
    package - package source - download
