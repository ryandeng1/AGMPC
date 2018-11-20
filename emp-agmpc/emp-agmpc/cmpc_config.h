#ifndef __CMPC_CONFIG
#define __CMPC_CONFIG
const static int abit_block_size = 1024;
const static int fpre_threads = 1;
#define LOCALHOST
//#define __MORE_FLUSH
//#define __debug
/*
const static char *IP[] = {""
,	"127.0.0.1"
,	"127.0.0.1"
,	"127.0.0.1"};
*/
const static char *IP[] = {"ec2-52-39-162-238.us-west-2.compute.amazonaws.com", 
"ec2-34-223-215-198.us-west-2.compute.amazonaws.com", 
"ec2-23-20-124-131.compute-1.amazonaws.com", 
"ec2-52-73-142-253.compute-1.amazonaws.com"};

const static bool lan_network = false;
#endif// __C2PC_CONFIG
