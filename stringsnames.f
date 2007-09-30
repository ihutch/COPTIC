c****************************************************************
c String construction, length, and concatenation routines.
c***************************************************************
C Construct a name based on a series of parameters.
c Example
c      program nameconstruct
c      character*30 name
cc Zero the name first. Very Important!
c      name=' '
c      call nameappendexp(name,'F',11.,2)
c      call nameappendint(name,'second',37,2)
c      np=nbcat(name,' dat')
c      write(*,*)name
c      end
c*******************************************************************
      subroutine nameappendexp(name,paramname,param,ip)
      character*(*) name,paramname
      real param
      character*12 sval
      sval=' '
      np=nbcat(name,paramname)      
      call eformati(param,sval,ip)
      np=nbcat(name,sval)
      end
c*******************************************************************
      subroutine nameappendint(name,paramname,iparam,ip)
      character*(*) name,paramname
      character*12 form
      integer iparam
      if(iparam.ge.0)then
         write(form,301)ip,ip
      else
         write(form,301)ip+1,ip
      endif
 301  format('(i',i1,'.',i1,')')
      np=nbcat(name,paramname)
      write(name(lentrim(name)+1:),form)iparam
      end
c******************************************************************
c Concatenate string2 trimmed of trailing blanks
c to string1 trimmed of trailing blanks. 
c Return the number of catenated characters.
      function nbcat(string1,string2)
      character*(*) string1,string2
      l1=lentrim(string1)
      l2=lentrim(string2)
      l2=min(l2,len(string1)-l1)
      string1(l1+1:l1+l2)=string2(1:l2)
      nbcat=l2
      end
c******************************************************************
c Obtain the length of a string omitting trailing blanks.
      function lentrim(string)
      character*(*) string
      do i=len(string),1,-1
         if(string(i:i).ne.' ') goto 101
      enddo
      i=0
 101  lentrim=i
      end
c******************************************************************
c Construct an exp-formatted value with ip significant figures.
      subroutine eformati(value,string,ip)
      character*(*) string
      character*6 form
      write(form,201)ip,ip
 201  format('(i',i1,'.',i1,')')
      if(value.lt.1.e-10*10.**ip)then
         iti=0
         it2=0
      elseif(value.gt.0.9e9*10.**ip)then
         iti=.99999*10.**ip
         it2=9
      else
         iti=nint(alog10(value)-0.49)-(ip-1)
         it2=nint(value/10.**iti)
      endif
      write(string(1:),form)it2
      if(iti.lt.0) then
         string(1+ip:1+ip)='m'
         iti=-iti
      else
         string(1+ip:1+ip)='e'
      endif
      write(string(2+ip:2+ip),'(i1.1)')iti
      end
c******************************************************************
c Return the position (not pointer) of match in string
c Trailing blanks are ignored in either string.
      function istrstr(string,match)
      integer istrstr
      character*(*) string,match
      ls=lentrim(string)
      lm=lentrim(match)
      istrstr=0
      do i=1,ls
         do j=1,lm
            if(string(i+j-1:i+j-1).ne.match(j:j)) goto 101
         enddo
c Here when matched.
         istrstr=i
         goto 102
 101  enddo
 102  return
      end
