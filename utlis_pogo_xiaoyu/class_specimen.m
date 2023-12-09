classdef class_specimen
   properties
       SampleNumber
       Material
       ct=0
       ct44 = 0
       ct55 = 0
       ct66 = 0
       cl
       alpha_l
       alpha_t
       rho1
       h1
       original
   end
%    properties (Dependent)
%       Modulus
%    end
   methods
       
       function specimen=class_specimen(SampleNumber, Material, ct, cl, alpha_l, alpha_t, rho1, h1)
           specimen.SampleNumber=SampleNumber;
           specimen.Material=Material;
           specimen.ct=ct;
           specimen.ct44=0;
           specimen.ct55=0;
           specimen.ct66=0;
           specimen.cl=cl;
           specimen.alpha_l=alpha_l;
           specimen.alpha_t=alpha_t;
           specimen.rho1=rho1;
           specimen.h1=h1;
       end
       
       % assign the orthotropic shear wave velocities
       function specimen = orthotropic_shear_vel(specimen, ct44, ct55, ct66)
           specimen.ct44 = ct44;
           specimen.ct55 = ct55;
           specimen.ct66 = ct66;
       end
       
       function specimen = change_h(specimen, new_h)
           if isempty(specimen.original)
               specimen.original=specimen.h1;
           end
           specimen.h1=new_h;
       end
       
       function specimen=recovery_h1(specimen)
           specimen.h1=specimen.original;
       end
       
       % delamination layer (air)
       function specimen=change_toAir(specimen)
           specimen.cl = 340;
           specimen.rho1 = 1.225;
       end

   end
      
end

