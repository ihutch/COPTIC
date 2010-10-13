      subroutine helpusage()
 301  format(a,i5,a,i5)
 302  format(a,f8.3)
      write(*,301)'Usage: ccpic [switches]'
      write(*,301)'Parameter switches.'
     $     //' Leave no gap before value. Defaults indicated [ddd'
      write(*,301)' -ni   set No of particles/node; zero => unset.    ['
     $     ,n_part
      write(*,301)' -rn   set reinjection number at each step.        ['
     $     ,ninjcomp
      write(*,302)' -ri   set rhoinfinity/node => reinjection number. ['
     $     ,ripernode
      write(*,302)' -rx   set Edge-potl relax rate: 0=>off, 1=>immed. ['
     $     ,crelax
      write(*,302)' -dt   set Timestep.              [',dt,
     $     ' -da   set Initial dt accel-factor[',bdt
      write(*,302)' -ds   set subcycle fraction.     [',subcycle
      write(*,301)' -s    set No of steps.           [',nsteps
      write(*,302)' -v    set Drift velocity.        [',vd
      write(*,302)' -t    set Ion Temperature.       [',Ti
      write(*,302)' -l    set Debye Length.          [',debyelen
      write(*,301)' -a    set averaging steps.       [',iavesteps
     $     ,'     Also period of diagnostic writes.'
      write(*,301)' -w    set write-step period.     [',iwstep
      write(*,301)' -m    set No of diag-moments(7). [',ndiags
      write(*,301)' -ct   set collision time.        [',colntime
      write(*,301)' -vn   set neutral drift velocity [',vneutral
c      write(*,301)' -xs<3reals>, -xe<3reals>  Set mesh start/end.'
      write(*,301)' -of<filename>  set name of object data file.'
     $     //'   [ccpicgeom.dat'
      write(*,301)
     $     ' -fs[path]  Attempt to restart from state saved [in path].'
      write(*,301)'Debugging switches for testing'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,301)' -gs[] Plot slices of solution potential, density. '
     $     //'[At step n]. [',ipstep
      write(*,301)' -gd -gp Turn off slicing of density, potential. '
      write(*,301)' -gf   set quantity plotted for flux evolution and'//
     $     ' final distribution. [',ifplot
      write(*,301)' -gc   set wireframe [& stencils(-)] mask.'//
     $     ' objects<->bits. [',iobpl
      write(*,301)' -gr   set wireframe override plot scale'//
     $     ' for -gc plot.  [',rcij
      write(*,301)' -go   set No of orbits'
     $     //'(to plot on objects set by -gc). [',norbits
      write(*,301)' -at   set test angle.'
     $     //' -an   set No of angles. '
      write(*,301)' -ck   set checking timestep No. [',ickst
      write(*,301)' -h -?   Print usage.'
      write(*,301)' -ho     Print geomobj file format description'
      call exit(0)
      end

