Troubleshooting and frequently asked questions (FAQs)
------------

Many errors that may be encountered may ultimately be the result of user error. If you encounter an error message any time that this pipeline is used, carefully check the command you used for any spelling errors. Additionally, many of these error messages give some detail as too where the code is wrong. Here are a few common errors and our suggestions for basic troubleshooting.

* Are you using the correct "profile" to run AmrPlusPlus?
  * We provide many examples of profile configurationg and choosing the correct one depends on your computing environment.
    * If you have singularity or installed on your server, we recommend using the "singularity" or profile to avoid the installation of any additional tools. 
    * If you already have the tools installed on your server, the best option is to configure the local.config file to point to the absolute PATH to each too.
* Are the right user permissions are granted to the file/directory/server in which you are going to run the pipeline?
  * In servers with multiple users, there are often cases in which certain directories give some users more editing privileges than others. Start by navigating to the directory in which you will be working. Next, type “ls -lha or ls -l”. This produces a list of all files in that directory and info on what permissions the user has using the “-rwxrwxrwx” scheme; r = read permissions, w = writing permissions, and x = execute permissions).
  * Permission errors could be due to the directories chosen for the pipeline output or individual bioinformatic tools installed by other users, for example. 
  * Review this tutorial for more information regarding file permissions: https://www.guru99.com/file-permissions.html

