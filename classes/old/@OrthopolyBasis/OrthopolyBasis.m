classdef OrthopolyBasis

  methods
    function obj = OrthopolyBasis(varargin)
    % OrthopolyBasis -- Class constructor for an Orthogonal Polynomial Basis
    %
    % obj = OrthopolyBasis({recurrence, normalization, weight_normalization, dim,
    %                       scale, shift, weight_function, description})
    %
    %     Creates a tensor-product orthogonal polynomial basis object. The
    %     optional inputs have the following default values:
    %
    %         recurrence:     The three-term recurrence relation. Is a function
    %                         handle with two outputs: [a,b] = recurrence(N),
    %                         and the outputs a and b are the first N
    %                         coefficients of the recurrence vectors:
    %                         p_{n+1} = (x - a_n) p_n - b_n p_{n-1},
    %                         and naturally a(1) = a_0, b(1) = b_0.
    %         normalization:  Specifies the normalization used to define the
    %                         polynomials. Standard values are:
    %                         'normal': the polynomials are L^2_w orthonormal,
    %                                   with w being the weight function.
    %                         'monic':  the polynomials are monic
    %         weight_normalization: The normalization used for the weight
    %                               function. The default is '' in which case
    %                               nothing is done. Other possibilities are
    %                               'probability', which means that the weight
    %                               function has unity L^1 norm. This affects,
    %                               among other things, the Gauss-type
    %                               quadrature rules and the 'normal'
    %                               polynomials.
    %        dim:             The tensor-product spatial dimension. Note that in
    %                         this case *all* of the optional arguments can be
    %                         dim x 1 cell arrays that contain the necessary
    %                         information for each dimension. If any of the
    %                         inputs is not a cell array, or is a cell array of
    %                         length 1, then that single input is repeated for
    %                         each dimension.
    %        scale:           The given recurrence relation is assumed to hold
    %                         for some standard interval I on which the
    %                         polynomials are defined. `scale' and `shift'
    %                         change this interval to scale*I + shift (for each
    %                         dimension). This does not affect the values of the
    %                         polynomials, but does affect quantities that rely
    %                         on the values of norms (such as the quadrature
    %                         rules).
    %        shift:           A corresponding affine parameter with meaning
    %                         explained in 'scale'.
    %        weight_function: A function handle returning the value of the
    %                         weight function on the standard interval I. This
    %                         is not necessary for any `standard' routines, but
    %                         some specialized ones do require knowledge of the
    %                         weight function, which cannot be divined from the
    %                         recurrence coefficients.
    %        description:     A struct of non-standard makeup that describes
    %                         what kind of polynomials these are ('chebyshev',
    %                         'hermite', etc.). Although this has non-standard
    %                         definition, some advanced or specialized routines
    %                         make use of certain fields in this struct.

      obj = parse_constructor_inputs(varargin{:});

  end 
end
