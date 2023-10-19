# source("./init_R.R")

# This script 
# - initializes a new renv project from r/4.2.1 or activate the already initialized one
# - installs packages from R_requirements.txt
# - reports status of the renv project
# - takes a snapshot of the renv project (to create a lockfile for reproducibility)

# initiate new renv project or restore existing one
if (!file.exists(".Rprofile")) {

    library(renv)

    # # disable automatic snapshotting
    # auto.snapshot <- getOption("renv.config.auto.snapshot")
    # options(renv.config.auto.snapshot = FALSE)

    renv::init(bare=TRUE)

    # install BiocManager beforehand
    renv::install("BiocManager@1.30.20")

    # read R packages from pseudo-requirements file: R_requirements.txt
    # and install packages. Line format is: package_name@version
    packages <- readLines("R_requirements.txt")
    for (package in packages) {
        # skip the lines starting with #
        if (substr(package, 1, 1) == "#") {
            next
        } 

        # try catch block to install packages
        tryCatch({
            renv::install(package)
        }, error = function(e) {
            print(paste("Searching ", package, " in bioconductor", sep=""))

            tryCatch({
                bioc_package <- strsplit(package, "@")[[1]][1]
                bioc_package <- paste("bioc::", bioc_package, sep="")
                renv::install(bioc_package)
                print(paste("Installed ", bioc_package, " from  bioconductor", sep=""))
            }, error = function(e) {
                print(paste("Could not find ", package, " in CRAN or bioconductor", sep=""))
                browser()
            })
        })
    }


    # take snapshot of the project
    renv::snapshot()


    # # restore automatic snapshotting
    # options(renv.config.auto.snapshot = auto.snapshot)


} else { # TODO add snapshot code to here

    library(renv)

    # list r packages
    renv::dependencies()

}

browser()



# take a look at the status
renv::status()

