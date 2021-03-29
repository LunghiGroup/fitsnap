        module meta_class
        use lists_class
        implicit none

         type memory_pot
          type(list)           :: gauss
          double precision     :: weight=1.0
          double precision     :: sigma=0.2
         end type memory_pot

         type(memory_pot)                :: meta1,meta2
         logical                         :: do_meta=.false.
         double precision                :: dist1,dist2,height1,height2
         double precision                :: f1(3),f2(3),f3(3)

        end module meta_class
