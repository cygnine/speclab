classdef PolynomialNormalization < Singleton
  properties
    description;
  end

 methods(Access=private)
    % Guard the constructor against external invocation.  We only want
    % to allow a single instance of this class.  See description in
    % Singleton superclass.
    function self = PolynomialNormalization()
       % Initialise your custom properties.
       self.description = 'polynomial normalization';
    end
 end
 
 methods(Static)
    % Concrete implementation.  See Singleton superclass.
    function self = instance()
       persistent PolynomialNormalization_instance
       if isempty(PolynomialNormalization_instance)
          self = PolynomialNormalization();
          PolynomialNormalization_instance= self;
       else
          self = PolynomialNormalization_instance;
       end
    end
 end
 
end
