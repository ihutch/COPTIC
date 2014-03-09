      program sysfonttest
c Uses native PS fonts in the file output
c while still vector fonts on screen.
      
      character*8 s1
      character*20 s2
      character*80 str
      character*2 str1

      call pfset(3)
c Set to native PS fonts.
      write(*,*) 'calling pfPSset'
      call pfPSset(1)
      write(*,*) 'calling pltinit'
      call pltinit(0.,1.,0.,1.)
      write(*,*) 'called pltinit'
c      call PSfontinit()
      call vecn(0.2,0.5,1)
      call drcstr('Blah abcdefg!dwxy!dz')
c      call charangl(-20.)
      call drcstr('Second !BdrcstrBCDS!@')
      call vecn(0.4,0.6,1)
      call charsize(.05,.05)
      call drcstr('Larger Font')

      call charsize(.03,.03)
      s2='ghijklmnopqrstuvwxyz'
      call vecn(.2,.3,0)
      call drcstr(s2)
      call vecn(.2,.25,0)
      call drcstr('!A'//s2)

      call charangl(50.)
      call drwstr(.2,.35,'50 degrees')

      call charangl(0.)

      call vecn(.5,.45,0)
      call jdrwstr(.6,.4,'Jdrwstr',0.)

      call vecn(0.6,0.4,1)

c      call gtic(2.43,1,.15,.15,1.,0.,1.,0.,.true.)
c      call axis()

      goto 1
c From here on is direct PS calls that won't plot on screen.
      call PSstartfont('/Helvetica')
      s1='abcdefgh'
      call PSchardrw(s1)

      call vecn(0.1,0.1,1)
      call PSchardrw(s2)
      call PSsetfont(1)
      call charsize(.03,.03)
      call charangl(60.)
      call PSchardrw(s2)
      call PSsetfont(2)
      call PSchardrw(s2)
      call PSsetfont(0)
      call PSchardrw(s2)

 
 1    call pltend()

c Now we try a case where we switch to PSfonts half way through

      call charsize(0.02,0.02)
      call pltinit(0.,1.,0.,1.)
      call axis()
      call pfPSset(1)
      call axlabels('Postscript','')
      call pfPSset(0)
      call axlabels('','Vector')
      call pfPSset(3)
      call jdrwstr(.5,.5,'Postscript, vector',0.)
      call pfPSset(3)
      call jdrwstr(.5,.45,'Postscript, vector',0.)
      call pltend()

      end
