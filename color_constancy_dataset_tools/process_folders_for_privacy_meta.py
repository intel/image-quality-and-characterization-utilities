#
# Copyright (c) 2019, Intel Corporation
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Intel Corporation nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import os, re
rootdir = '/example_path'

existing_meta_files = 0
new_meta_files_created = 0

# Iterate through the folders in the 'rootdir' parent folder
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
	if ".plain16" in file:
		# Update the .meta file
		meta_file = os.path.join(subdir, file)
		meta_file = meta_file.replace('.plain16','.meta')
		
		exists = os.path.isfile(meta_file)
		if exists:
			out_str = ""
			tag_found = 0
			with open(meta_file,'r+') as file2:
				for line in file2:
					if 'privacy_masked' in line:
						out_str = out_str + "privacy_masked\t1\n"
						tag_found = 1
					elif not line.isspace():
						out_str = out_str + line

			if tag_found < 1:
				out_str = out_str + "privacy_masked\t1\n"
				f = open(meta_file,"w+")
				f.write(out_str)
				f.close()
			existing_meta_files = existing_meta_files + 1
			print ".meta file already exists"
		else:
			out_str = "privacy_masked\t1\n"
			f = open(meta_file,"w+")
			f.write(out_str)
			f.close()
			new_meta_files_created = new_meta_files_created + 1
			print "Created new .meta file with: privacy_masked	1"		

# Print statistics of the processed files
print "Existing .meta files: " + str(existing_meta_files)
print "New .meta files created: " + str(new_meta_files_created)

