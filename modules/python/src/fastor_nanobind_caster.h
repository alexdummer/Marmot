#pragma once

#include <Fastor/Fastor.h>
#include <algorithm>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace nanobind::detail {
  template < typename T, size_t... Dims >
  struct type_caster< Fastor::Tensor< T, Dims... > > {
    using TensorType = Fastor::Tensor< T, Dims... >;
    NB_TYPE_CASTER( TensorType, const_name( "numpy.ndarray[" ) + make_caster< T >::Name + const_name( "]" ) );

    bool from_python( handle src, uint8_t flags, cleanup_list* cleanup ) noexcept
    {
      using NdArray = nb::ndarray< T, nb::shape< Dims... >, nb::c_contig, nb::device::cpu >;

      type_caster< NdArray > nd_caster;
      if ( !nd_caster.from_python( src, flags, cleanup ) ) {
        return false;
      }

      NdArray nd = nd_caster.value;
      std::copy( nd.data(), nd.data() + nd.size(), value.data() );
      return true;
    }

    static handle from_cpp( const TensorType& src, rv_policy policy, cleanup_list* cleanup ) noexcept
    {
      size_t shape[] = { Dims... };
      size_t ndim    = sizeof...( Dims );
      size_t size    = 1;
      for ( size_t i = 0; i < ndim; ++i )
        size *= shape[i];

      T* data = new T[size];
      std::copy( src.data(), src.data() + size, data );

      nb::capsule owner( data, []( void* p ) noexcept { delete[] (T*)p; } );

      // Use generic nanobind cast to create python object from ndarray
      nb::object obj = nb::cast( nb::ndarray< nb::numpy, T >( data, ndim, shape, owner ) );
      return obj.release();
    }
  };
} // namespace nanobind::detail
