# misc
# show unique string lengths in a table column (here: column 7 of file mappingtests.tab)
cut -f7 mappingtests.tab | awk '{print length}' | sort | uniq

# du -> (d)isk (u)sage
# sizes of folders and subfolders
du -sh /path/to/directory		# (s)ummarize subfolders in indicated folder / (h)uman readable
du /path/to/directory | sort -n -r	# sort by size in (r)everse (n)umerical order
du --max-depth=2 /path/to/directory	# summarize subfolders more than two levels below indicated folder

# df -> (d)isk (f)ree
# free space on available file systems
df -h					# (h)uman readable

# tar -> (t)ape (ar)chive
tar -xvf <file>		
tar -xvzf <file>	# tar.gz

# cp files from subdirectories of source folder to destination folder
#	* = file filter
#	/Y = suppress overwrite confirmation
for /R "[source folder]" %f in (*) do copy /Y "%f" "[destination folder]"

# create .tar.gz
tar -pczf name_of_your_archive.tar.gz /path/to/directory	# (p)ermissions, (f)ile, (c)reate, g(z)ip
# untar .tar.gz
tar -xzvf name_of_archive /extract/here				# e(x)tract, g(z)ip, (v)erbose, (f)ile
