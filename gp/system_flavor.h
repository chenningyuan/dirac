#define LINUX   (0)
#define FREEBSD (1)
#define OPENBSD (2)
#define NETBSD  (3)
#define DARWIN  (4)
#define UNKNOWN (-1)

#if ! defined SYS_FLAVOR
#define SYS_FLAVOR UNKNOWN 
#endif
