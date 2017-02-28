PATH ORGANISATION
    bin : version specific executables *
    doc : version specific documentation * 
    log : experiment specific log files  **
    raw_data : experiment specific raw_data (currently --> Workflow input) 
    results : experiment specific results **
    run : all files related to this version specific workflow (currently --> mapping id files) 
    src : version specific scripts 
    tmp : experiment specific temporary files **
    
    
* commited with .gitkeep ==> for every commit, commit this folder regardless of its empty or not.
** commited with .gitignore ==> for every commit, ignore files contained in this folder.

Folder not flagged with "*" or "**" do not contain .gitignore or .gitkeep files because they are not empty or/and dont needed it.
