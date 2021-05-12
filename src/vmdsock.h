/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2008 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: vmdsock.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $      $Date: 2008/03/27 19:36:52 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   socket interface layer, abstracts platform-dependent routines/APIs
 ***************************************************************************/

#if defined(VMDSOCKINTERNAL)

#if !defined(_MSC_VER)
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <sys/file.h>
#endif

typedef struct {
  struct sockaddr_in addr; /* address of socket provided by bind() */
  int addrlen;             /* size of the addr struct */
  int sd;                  /* socket file descriptor */
} vmdsocket;

#endif /* VMDSOCKINTERNAL */

#include "molmodel/internal/Linkage.h"

#ifdef __cplusplus
extern "C" {
#endif

SimTK_MOLMODEL_EXPORT  int  vmdsock_init(void);
SimTK_MOLMODEL_EXPORT  void *vmdsock_create(void);
SimTK_MOLMODEL_EXPORT  int  vmdsock_bind(void *, int);
SimTK_MOLMODEL_EXPORT  int  vmdsock_listen(void *);
SimTK_MOLMODEL_EXPORT  void *vmdsock_accept(void *);  /* return new socket */
SimTK_MOLMODEL_EXPORT  int  vmdsock_connect(void *, const char *, int);
SimTK_MOLMODEL_EXPORT  int  vmdsock_write(void *, const void *, int);
SimTK_MOLMODEL_EXPORT  int  vmdsock_read(void *, void *, int);
SimTK_MOLMODEL_EXPORT  int  vmdsock_selread(void *, int);
SimTK_MOLMODEL_EXPORT  int  vmdsock_selwrite(void *, int);
SimTK_MOLMODEL_EXPORT  void vmdsock_shutdown(void *);
SimTK_MOLMODEL_EXPORT  void vmdsock_destroy(void *);

#ifdef __cplusplus
}
#endif

