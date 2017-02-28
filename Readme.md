GENIAL : gene set enrichment analysis
===================


This is a quick description of GENIAL tools. Don't delete me, I'm very helpful! I can be recovered anyway in the **Utils** tab of the <i class="icon-cog"></i> **Settings** dialog.

----------




<i class="icon-file"></i>  Path organisation
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

