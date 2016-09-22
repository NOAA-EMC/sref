      subroutine uvoedz(ulnod,vlnev,flnod,flnev,epse,epso,
     x                  ls_node)
cc
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_physcons, rerth => con_rerth
      implicit none
cc
      real(kind=kind_evod) ulnod(len_trio_ls,2)
      real(kind=kind_evod) vlnev(len_trie_ls,2)
      real(kind=kind_evod) flnod(len_trio_ls,2)
      real(kind=kind_evod) flnev(len_trie_ls,2)
cc
      real(kind=kind_evod)  epse(len_trie_ls)
      real(kind=kind_evod)  epso(len_trio_ls)
cc
      integer              ls_node(ls_dim,3)
cc
cc
      integer              l,locl,n
cc
      integer              indev,indev1,indev2
      integer              indod,indod1,indod2
      integer              inddif
cc
      real(kind=kind_evod) rl,rn,rnp1
cc
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons2     !constant
cc
      integer              indlsev,jbasev
      integer              indlsod,jbasod
cc
      include 'function2'
cc
cc......................................................................
cc
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
cc
cc
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L) - 1
         else
            indev2 = indlsev(jcap  ,L) - 1
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
cc
         rl   = l
         rn   = l+1
         rnp1 = l+1+1
         do indev = indev1 , indev2
cc
                                          flnod(indev-inddif,1) =
     x                              -rl * ulnod(indev-inddif,2)
     x      + rn   * epse(indev+1     ) * vlnev(indev+1     ,1)
     x      - rnp1 * epso(indev-inddif) * vlnev(indev       ,1)
cc
                                          flnod(indev-inddif,2) =
     x                               rl * ulnod(indev-inddif,1)
     x      + rn   * epse(indev+1     ) * vlnev(indev+1     ,2)
     x      - rnp1 * epso(indev-inddif) * vlnev(indev       ,2)
cc
              rn   = rn   + cons2     !constant
              rnp1 = rnp1 + cons2     !constant
         end do
cc
      end do
cc
cc......................................................................
cc
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
cc
           rl = l
           rn = l
         do indev = indev1 , indev2
cc
                                        flnev(indev       ,1) =
     x                            -rl * vlnev(indev       ,2)
     x      - rn * epso(indev-inddif) * ulnod(indev-inddif,1)
cc
                                        flnev(indev       ,2) =
     x                             rl * vlnev(indev       ,1)
     x      - rn * epso(indev-inddif) * ulnod(indev-inddif,2)
cc
              rn = rn + cons2     !constant
         end do
cc
      end do
cc
cc......................................................................
cc
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(l,l) + 1
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap-1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         indod1 = indlsod(l+1,l)
         inddif = indev1 - indod1
cc
         rnp1 = l+2+1
         do indev = indev1 , indev2
cc
                                   flnev(indev       ,1) =
     x                             flnev(indev       ,1)
     x      + rnp1 * epse(indev) * ulnod(indev-inddif,1)
cc
                                   flnev(indev       ,2) =
     x                             flnev(indev       ,2)
     x      + rnp1 * epse(indev) * ulnod(indev-inddif,2)
cc
              rnp1 = rnp1 + cons2     !constant
         end do
cc
      end do
cc
cc......................................................................
cc
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
cc
         if (mod(L,2).eq.mod(jcap+1,2)) then
cc          set the even (n-l) terms of the top row to zero
            flnev(indlsev(jcap+1,l),1) = cons0     !constant
            flnev(indlsev(jcap+1,l),2) = cons0     !constant
         else
cc          set the  odd (n-l) terms of the top row to zero
            flnod(indlsod(jcap+1,l),1) = cons0     !constant
            flnod(indlsod(jcap+1,l),2) = cons0     !constant
         endif
cc
cc
      enddo
cc
cc......................................................................
cc
cc
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         jbasod=ls_node(locl,3)
         indev1 = indlsev(L,L)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
            indod2 = indlsod(jcap  ,L)
         else
            indev2 = indlsev(jcap  ,L)
            indod2 = indlsod(jcap+1,L)
         endif
         do indev = indev1 , indev2
            flnev(indev,1)=flnev(indev,1)/rerth
            flnev(indev,2)=flnev(indev,2)/rerth
         enddo
cc
         do indod = indod1 , indod2
            flnod(indod,1)=flnod(indod,1)/rerth
            flnod(indod,2)=flnod(indod,2)/rerth
         enddo
cc
cc
      enddo
cc
      return
      end
