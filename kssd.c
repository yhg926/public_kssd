//	Copyright 2019 Huiguang Yi. All Rights Reservered.
//
//	Licensed under the Apache License, Version 2.0 (the "License");
//	you may not use this file except in compliance with the License.
//	You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//	Unless required by applicable law or agreed to in writing, software
//	distributed under the License is distributed on an "AS IS" BASIS,
//	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//	See the License for the specific language governing permissions and
//	limitations under the License.
   
#include "global_wrapper.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_basic.h"
const char *argp_program_version = "KSSD version 1.2.21";
const char *argp_program_bug_address = "yhg926@gmail.com";
int main(int argc, char** argv)
{
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);
  setvbuf(stderr, (char *)NULL, _IOLBF, 0);
 const char subcommand_n[] = "<subcommand>";
  domain = (char*) malloc(strlen(argv[0])+1);
  strcpy(domain,argv[0]);
  long_domain = (char*) malloc( strlen(domain) + strlen(subcommand_n) + 2);
  snprintf(long_domain,strlen(domain) + strlen(subcommand_n) + 2,"%s %s",domain,subcommand_n);
  cmd_global(argc, argv);
  return 0;
}
