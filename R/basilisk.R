#' @importFrom basilisk BasiliskEnvironment
env_macs2 <- BasiliskEnvironment("env_macs2", pkgname="peakPeekeR",
                                 packages=c("macs2==2.2.7.1", "python==3.8"),
                                 channels=c("bioconda", "conda-forge"))

env_macs <- BasiliskEnvironment("env_macs", pkgname="peakPeekeR",
                                packages=c("macs==1.4.3", "python==2.7"),
                                channels=c("bioconda", "conda-forge"))

env_sicer2 <- BasiliskEnvironment("env_sicer2", pkgname="peakPeekeR",
                                packages=c("sicer2==1.0.3", "python==3.8"),
                                channels=c("bioconda", "conda-forge"))

env_genrich <- BasiliskEnvironment("env_genrich", pkgname="peakPeekeR",
                                  packages=c("genrich==0.6.1", "python==3.8"),
                                  channels=c("bioconda", "conda-forge"))
