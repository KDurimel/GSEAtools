GENIAL : gene set enrichment analysis
===================


This is a quick description of GENIAL tools path organization. Don't delete me, I'm very helpful!

For more details about how GENIAL tools work and how to use it, please go to the documentation folder (doc).

----------




<i class="icon-file"></i>  Path organization
=================

    bin : version specific executables [1]
    doc : version specific documentation[1]
    log : experiment specific log files[2]
    **raw_data : experiment specific raw_data (currently --> Workflow input) 
    **results : experiment specific results**[2]
    run : all files related to this version specific workflow 
    (currently --> mapping   id files) 
    src : version specific scripts 
    **tmp : experiment specific temporary files**[2]


 


----------


  [1]: Files commited with [gitkeep](http://www.murderdev.net/2012/10/gitignore-gitkeep/) : for every commit, ignore files contained in this folder.
   
   [2]: Files commited with [gitignore](https://git-scm.com/docs/gitignore) : for every commit, commit this folder regardless of its empty or not.

Folders who are not flagged with "*" or "**" do not contain .gitignore or .gitkeep files because they are not empty or/and dont needed it.

