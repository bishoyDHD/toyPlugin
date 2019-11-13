# trekC
minimal cooker: Analysis framework/toolkit compatible with ROOT v6 courtesy of Jan C. Bernauer

FULL DISCLAIMER:
----------------
Installing mini-Cooker is a lot like hanging out with a pubescent child, ANNOYING! Once installed however, there is a lot of functionality and usefulness limited only by the user's imagination.

HOW TO INSTALL:
---------------
Prerequisites:
* Need ~~[GenFit](https://github.com/GenFit/GenFit) and [eigen3](https://github.com/eigenteam/eigen-git-mirror)~~
* Usual: GSL, xqilla, lapack, expat
* Boost between 1.47-1.62 (will conflict with std::cpp if higher)
* CLHEP (might be system dependent feature, see CLHEP section below)
* Need C++-11 (fortunately this is set in CMakeLists.txt)

### Installing CLHEP:
With respect to my machine I found that my install dir CLHEP could not be /opt/clhep (recommended installation dir)-- otherwise the include libs cannot be accessed/found (even when the environment variables are set). For my case this makes sense because I am running on rather secure OS (same case with bleeding-edge Os's). The way to overcome this difficulty/inconvenience was to *cmake* without explicitly setting the **install dir**, which means that the build dir and **install dir** are the same. The benefit is that the **include dir** become /usr/include, which solves the problem. Again, this depends on your OS of course.
