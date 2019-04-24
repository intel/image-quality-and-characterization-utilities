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

deleted_plain16_cnt = 0
renamed_privacy_plain16_cnt = 0
updated_meta_cnt = 0
new_meta_cnt = 0
deleted_privacy_jpg_cnt = 0

# Iterate through the folders in the 'rootdir' parent folder
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
		if "_privacy.plain16" in file:
			# Replace the .plain16 raw file with the privacy masked version
		    	print "Modifying: " + os.path.join(subdir, file)
			deleted_file = file.replace('_privacy','')
			print "Deleted file: " + os.path.join(subdir, deleted_file)
			try:
				os.remove(os.path.join(subdir, deleted_file))
				deleted_plain16_cnt = deleted_plain16_cnt + 1
			except:
				print "Could not delete file " + os.path.join(subdir, deleted_file)
			os.rename(os.path.join(subdir, file), os.path.join(subdir, deleted_file))
			renamed_privacy_plain16_cnt = renamed_privacy_plain16_cnt + 1

			# Update the .meta file
			meta_file = os.path.join(subdir, deleted_file)
			meta_file = meta_file.replace('.plain16','.meta')
			out_str = ""
			tag_found = 0
			try:
				with open(meta_file,'r+') as file2:
					for line in file2:
						if 'privacy_masked' in line:
							out_str = out_str + "privacy_masked\t1\n"
							tag_found = 1
						elif not line.isspace():
							out_str = out_str + line
				updated_meta_cnt = updated_meta_cnt + 1
			except:
				print "No .meta file found"
				new_meta_cnt = new_meta_cnt + 1
			if tag_found < 1:
				out_str = out_str + "privacy_masked\t1\n"
			f = open(meta_file,"w+")
			f.write(out_str)
			f.close()
		
			# Delete the _privacy.jpg file if there is one
			deleted_jpgfile = file.replace('.plain16','.jpg')		
			try:
				os.remove(os.path.join(subdir, deleted_jpgfile))
				print "Deleted jpg file: " + os.path.join(subdir, deleted_jpgfile)
				deleted_privacy_jpg_cnt = deleted_privacy_jpg_cnt + 1
			except:
				print "No _privacy.jpg file found"

# Find potential remaining _privacy.jpg files and remove them
for subdir, dirs, files in os.walk(rootdir):
	for file in files:
		if "_privacy.jpg" in file:
			os.remove(os.path.join(subdir, file))
			deleted_privacy_jpg_cnt = deleted_privacy_jpg_cnt + 1
			print "Deleted jpg file: " + os.path.join(subdir, file)

# Print statistics of the processed files
print "Deleted _plain16 count: " + str(deleted_plain16_cnt)
print "Renamed _privacy_plain16 count: " + str(renamed_privacy_plain16_cnt)
print "Updated meta files count: " + str(updated_meta_cnt)
print "New metas count: " + str(new_meta_cnt)
print "Deleted _privacy.jpg count: " + str(deleted_privacy_jpg_cnt)

