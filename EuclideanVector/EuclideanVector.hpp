//
//   _____               _   _              _                  _____
//  |  ___|             | | | |            | |                /  __ \   _      _
//  | |__   _   _   ___ | | | |  ___   ___ | |_   ___   _ __  | /  \/ _| |_  _| |_
//  |  __| | | | | / __|| | | | / _ \ / __|| __| / _ \ | '__| | |    |_   _||_   _|
//  | |___ | |_| || (__ \ \_/ /|  __/| (__ | |_ | (_) || |    | \__/\  |_|    |_|
//  \____/  \__,_| \___| \___/  \___| \___| \__| \___/ |_|     \____/
//
// 
//	Git Hub : https://github.com/Rasukairosu/EuclideanVector.git
//	
//	2024/10/22~
//
// 	Author : Hiyu Taninaka
// 
//	Ver 2.0	(Last Update 2024/10/31)
// 
//	-Japanese-
// 
//	EuclideanVector」シリーズは、1～4次元のベクトルを表現するクラスです。
//	計算結果の出力に集成体を使用することで、高速な計算を可能にしています。
//
//	基本的に組み込み算術型を前提としていますが、任意の型をメンバ型として指定することも可能です。
//	ただし、指定する型には引数なしのコンストラクタが必要で、const修飾子や参照型は使用できません。
//	また、演算子のオーバーロードについても、指定型で該当演算子が呼び出し可能な場合にのみ実体化されます。
//
//	メンバ関数も指定型での演算子の利用可否に依存して使用可能かが決まるため、この点には注意が必要です。
//	演算子の結果として返される型は EuclideanVector ではなく、計算専用の ResultPack 型であることにも留意してください。
//	この型はメンバ変数に直接アクセスできます。
//
//	EuclideanVector には EuclideanRecVector と EuclideanCmplVector の二種類があります。
//
//	1 : EuclideanRecVector
// 
//		データを再帰的に表現するため、自身より低次元のベクトルへの変換が可能です。
//		三次元ベクトルではxy()、四次元ベクトルではxyz()を使って変換してください。
//		構築に時間がかかるため、例えばfor文内でのローカル変数としての利用は推奨しません。
//		メンバ変数への直接アクセスはできないため、アクセスには x() や y() といったアクセサを使用してください。
//		値を一括で設定する場合は、v = EuclideanVector<T>(...) よりも set(...) を使用することでオーバーヘッドを削減できます。
//
//	2 : EuclideanCmplVector
// 
//		全要素を単独で保持する構造であるため、低次元のベクトルへの変換はできませんが、構築が高速です。
//		また、メンバ変数への直接アクセスも可能です。
// 
//	-English-
//
//	The "EuclideanVector" series represents vectors in 1 to 4 dimensions,
//	optimized for high - speed calculations by using aggregates for result output.
//
//	While built - in arithmetic types are the default, 
//	custom types can also be specified as member types.
//	These specified types must have a default constructor and cannot be const or reference types.
//	Operators are instantiated only if the specified type supports the corresponding operator.
//
//	Similarly, member functions are conditionally available,
//	depending on the specified type’s operator support.
//	Note that the result of an operation is not returned as a EuclideanVector type,
//	but rather as a ResultPack type, which provides direct access to its member variables.
//
//	The EuclideanVector series includes two types : EuclideanRecVector and EuclideanCmplVector.
//
//	1 : EuclideanRecVector
// 
//		This type uses a recursive representation of data, 
//		allowing conversion to lower - dimensional vectors.
//		For instance, use xy() for 3 - dimensional vectors and xyz() for 4 - dimensional vectors.
//		Due to its slower construction, 
//		using EuclideanRecVector as a local variable in a loop is not recommended.
//		Direct access to member variables is unavailable; instead, use accessor functions like x() and y().
//		For setting values all at once, use set(...) instead of v = EuclideanVector<T>(...) to reduce overhead.
//
//	2 : EuclideanCmplVector
// 
//		Unlike EuclideanRecVector, EuclideanCmplVector stores all elements independently.
//		While it does not support conversion to lower - dimensional vectors,
//		its construction is faster.It also allows direct access to member variables.
// 
//.

#ifndef THL_EUCLID_VECTOR_HPP
#define THL_EUCLID_VECTOR_HPP

#include <cmath>
#include <utility>
#include <type_traits>

#if defined(__clang__) || defined(__GNUC__)
#	define EUCVECTORINLINE inline
#	if (__cplusplus >= 	202002L)
#	define EUCNODISCARD				[[nodiscard]]
#	define EUCNODISCARD_MSG(msg)	[[nodiscard(msg)]]
#	elif (__cplusplus >= 201703L)
#	define EUCNODISCARD				[[nodiscard]]
#	define EUCNODISCARD_MSG(msg)	[[nodiscard]]
#	else
#	define EUCNODISCARD				
#	define EUCNODISCARD_MSG(msg)	
#	endif
#elif defined(_MSC_VER)
#	if (_MSC_VER  >= 1200)
#	define EUCVECTORINLINE __forceinline
#	else
#	define EUCVECTORINLINE inline
#	endif
#	if (_MSVC_LANG >= 	202002L)
#	define EUCNODISCARD				[[nodiscard]]
#	define EUCNODISCARD_MSG(msg)	[[nodiscard(msg)]]
#	elif (_MSVC_LANG >= 201703L)
#	define EUCNODISCARD				[[nodiscard]]
#	define EUCNODISCARD_MSG(msg)	[[nodiscard]]
#	else
#	define EUCNODISCARD				
#	define EUCNODISCARD_MSG(msg)	
#	endif
#endif

//name space begin.
namespace thl {
namespace vector {

//meta functions.
namespace meta {

	struct evd_euc_vec {};

	template<bool IF>
	using if_t = _STD enable_if_t<IF, int>;

	template<class T>
	using r_ref = _STD add_rvalue_reference_t<T>;

	template<bool IF>
	constexpr _STD add_rvalue_reference_t<meta::if_t<IF>> cnd_type() noexcept;

	template<class B, class T>
	using is_type_t = _STD enable_if_t<_STD is_same_v<B,B>, T>;

	template<class T, class B, bool IF>
	using co_inst_if_t = _STD enable_if_t<_STD is_same_v<T, B>&& IF, int>;

	template<class T>
	using no_ref = _STD remove_reference_t<T>;

	template<class C, class Head, class... Any>
	struct is_constructible_anynum_param 
		: private is_constructible_anynum_param<C, Any...>{
		static constexpr bool value = _STD is_constructible_v<C, Head> &&
			is_constructible_anynum_param<C, Any...>::value;
	};

	template<class C, class Head>
	struct is_constructible_anynum_param<C, Head> {
		static constexpr bool value = _STD is_constructible_v<C, Head>;
	};

	template<class Head, class... Body>
	struct num_of_templates 
		: private num_of_templates<Body...> {
		static constexpr size_t value = 1 + num_of_templates<Body...>::value;
	};

	template<class Head>
	struct num_of_templates<Head> {
		static constexpr size_t value = 1;
	};

	class not_equal_info {

		template<class L, class R>
		static constexpr auto can(int) -> decltype(_STD declval<L>() != _STD declval<R>(), _STD true_type());

		template<class L, class R>
		static constexpr _STD false_type can(...);

	public:
		template<class L, class R>
		static constexpr bool value = decltype(can<L, R>(0))::value;
	};

	class add_equal_info {

		template<class L, class R>
		static constexpr auto can(int) -> decltype(_STD declval<L>() += _STD declval<R>(), _STD true_type());

		template<class L, class R>
		static constexpr _STD false_type can(...);

	public:
		template<class L, class R>
		static constexpr bool value = decltype(can<L, R>(0))::value;
	};

	class sub_equal_info {

		template<class L, class R>
		static constexpr auto can(int) -> decltype(_STD declval<L>() -= _STD declval<R>(), _STD true_type());

		template<class L, class R>
		static constexpr _STD false_type can(...);

	public:
		template<class L, class R>
		static constexpr bool value = decltype(can<L, R>(0))::value;
	};

	class mul_equal_info {

		template<class L, class R>
		static constexpr auto can(int) -> decltype(_STD declval<L>() *= _STD declval<R>(), _STD true_type());

		template<class L, class R>
		static constexpr _STD false_type can(...);

	public:
		template<class L, class R>
		static constexpr bool value = decltype(can<L, R>(0))::value;
	};

	class div_equal_info {

		template<class L, class R>
		static constexpr auto can(int) -> decltype(_STD declval<L>() /= _STD declval<R>(), _STD true_type());

		template<class L, class R>
		static constexpr _STD false_type can(...);

	public:
		template<class L, class R>
		static constexpr bool value = decltype(can<L, R>(0))::value;
	};

	template<class To, class... Any>
	constexpr bool is_constructible_anynum_param_v = is_constructible_anynum_param<To, Any...>::value;

	template<class... Any>
	constexpr size_t num_of_templates_v = num_of_templates<Any...>::value;

	// !=.
	template<class Left, class Right>
	constexpr bool is_not_equal_v = not_equal_info::value<Left, Right>;

	// +=.
	template<class Left, class Right>
	constexpr bool is_add_equal_v = add_equal_info::value<Left, Right>;

	// -=.
	template<class Left, class Right>
	constexpr bool is_sub_equal_v = sub_equal_info::value<Left, Right>;

	// *=.
	template<class Left, class Right>
	constexpr bool is_mul_equal_v = mul_equal_info::value<Left, Right>;

	// /=.
	template<class Left, class Right>
	constexpr bool is_div_equal_v = div_equal_info::value<Left, Right>;

}

//details.
namespace detail {

	template<class E>
	struct ResultPacker_1 {
		E x;

		template<class T>
		EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
			EUCVECTORINLINE auto operator+(ResultPacker_1<T>&& pack) const noexcept(noexcept(ResultPacker_1<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}))
			->decltype(ResultPacker_1<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}) {
			return { _STD move(x) + _STD move(pack.x) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
			EUCVECTORINLINE auto operator-(ResultPacker_1<T>&& pack) const noexcept(noexcept(ResultPacker_1<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}))
			->decltype(ResultPacker_1<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}) {
			return { _STD move(x) - _STD move(pack.x) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE auto operator*(T&& scl) noexcept(noexcept(ResultPacker_1<meta::no_ref<decltype(x * scl)>>{x * scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_1<meta::no_ref<decltype(x* scl)>>{x* scl}) {
			return { x * scl };
		}

		template<class T>
		friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE static auto operator*(T&& scl, ResultPacker_1<E>&& right) noexcept(noexcept(ResultPacker_1<meta::no_ref<decltype(right.x * scl)>>{right.x * scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_1<meta::no_ref<decltype(right.x * scl)>>{right.x * scl}) {
			return { right.x * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
			EUCVECTORINLINE auto operator/(T&& scl) noexcept(noexcept(ResultPacker_1<meta::no_ref<decltype(x / scl)>>{x / scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_1<meta::no_ref<decltype(x / scl)>>{x / scl}) {
			return { x / scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator==(ResultPacker_1<T>&& right) noexcept(noexcept(bool(_STD move(x) == _STD move(right.x))))
			-> decltype(bool(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return _STD move(x) == _STD move(right.x);
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_1<T>&& right) noexcept(noexcept(bool(_STD move(x) != _STD move(right.x))))
			-> decltype(bool(_STD move(x) != _STD move(right.x)), _STD declval<bool>()) {
			return _STD move(x) != _STD move(right.x);
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_1<T>&& right) const noexcept(noexcept(bool(!(_STD move(x) == _STD move(right.x)))))
			-> decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), !(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return  !(_STD move(x) == _STD move(right.x));
		}

		void operator=(ResultPacker_1<E>) = delete;

	};

	template<class E>
	struct ResultPacker_2 {
		E x,y;

		template<class T>
		EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
			EUCVECTORINLINE auto operator+(ResultPacker_2<T>&& pack) const noexcept(noexcept(ResultPacker_2<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}))
			->decltype(ResultPacker_2<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}) {
			return { _STD move(x) + _STD move(pack.x), _STD move(y) + _STD move(pack.y) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
			EUCVECTORINLINE auto operator-(ResultPacker_2<T>&& pack) const noexcept(noexcept(ResultPacker_2<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}))
			->decltype(ResultPacker_2<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}) {
			return { _STD move(x) - _STD move(pack.x), _STD move(y) - _STD move(pack.y) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE auto operator*(T&& scl) noexcept(noexcept(ResultPacker_2<meta::no_ref<decltype(x * scl)>>{x * scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_2<meta::no_ref<decltype(x * scl)>>{x * scl}) {
			return { x * scl,y * scl };
		}

		template<class T>
		friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE static auto operator*(T&& scl, ResultPacker_2<E>&& right) noexcept(noexcept(ResultPacker_2<meta::no_ref<decltype(right.x* scl)>>{right.x* scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_2<meta::no_ref<decltype(right.x* scl)>>{right.x* scl}) {
			return { right.x * scl,right.y * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
			EUCVECTORINLINE auto operator/(T&& scl) noexcept(noexcept(ResultPacker_2<meta::no_ref<decltype(x / scl)>>{x / scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_2<meta::no_ref<decltype(x / scl)>>{x / scl}) {
			return { x / scl,y / scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator==(ResultPacker_2<T>&& right) noexcept(noexcept(bool(_STD move(x) == _STD move(right.x))))
			-> decltype(bool(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) == _STD move(right.x)) && (_STD move(y) == _STD move(right.y));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_2<T>&& right) noexcept(noexcept(bool(_STD move(x) != _STD move(right.x))))
			-> decltype(bool(_STD move(x) != _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) != _STD move(right.x)) || (_STD move(y) != _STD move(right.y));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_2<T>&& right) const noexcept(noexcept(bool(!(_STD move(x) == _STD move(right.x)))))
			-> decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), !(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return  (!(_STD move(x) == _STD move(right.x))) || (!(_STD move(y) == _STD move(right.y)));
		}

		void operator=(ResultPacker_2<E>) = delete;

	};

	template<class E>
	struct ResultPacker_3 {
		E x,y,z;

		template<class T>
		EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
			EUCVECTORINLINE auto operator+(ResultPacker_3<T>&& pack) const noexcept(noexcept(ResultPacker_3<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}))
			->decltype(ResultPacker_3<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}) {
			return { _STD move(x) + _STD move(pack.x), _STD move(y) + _STD move(pack.y), _STD move(z) + _STD move(pack.z) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
			EUCVECTORINLINE auto operator-(ResultPacker_3<T>&& pack) const noexcept(noexcept(ResultPacker_3<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}))
			->decltype(ResultPacker_3<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}) {
			return { _STD move(x) - _STD move(pack.x), _STD move(y) - _STD move(pack.y), _STD move(z) - _STD move(pack.z) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE auto operator*(T&& scl) noexcept(noexcept(ResultPacker_3<meta::no_ref<decltype(x * _STD forward<T>(scl))>>{x * _STD forward<T>(scl)}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_3<meta::no_ref<decltype(x * _STD forward<T>(scl))>>{x * _STD forward<T>(scl)}) {
			return { x * scl,y * scl,z * scl };
		}

		template<class T>
		friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE static auto operator*(T&& scl, ResultPacker_3<E>&& right) noexcept(noexcept(ResultPacker_3<meta::no_ref<decltype(right.x* scl)>>{right.x* scl}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_3<meta::no_ref<decltype(right.x* scl)>>{right.x* scl}) {
			return { right.x * scl,right.y * scl,right.z * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
			EUCVECTORINLINE auto operator/(T&& scl) noexcept(noexcept(ResultPacker_3<meta::no_ref<decltype(x / _STD forward<T>(scl))>>{x / _STD forward<T>(scl)}))
			->decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_3<meta::no_ref<decltype(x / _STD forward<T>(scl))>>{x / _STD forward<T>(scl)}) {
			return { x * scl,y * scl,z * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator==(ResultPacker_3<T>&& right) noexcept(noexcept(bool(_STD move(x) == _STD move(right.x))))
			-> decltype(bool(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) == _STD move(right.x)) && (_STD move(y) == _STD move(right.y)) && (_STD move(z) == _STD move(right.z));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_3<T>&& right) noexcept(noexcept(bool(_STD move(x) != _STD move(right.x))))
			-> decltype(bool(_STD move(x) != _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) != _STD move(right.x)) || (_STD move(y) != _STD move(right.y)) || (_STD move(z) != _STD move(right.z));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_3<T>&& right) const noexcept(noexcept(bool(!(_STD move(x) == _STD move(right.x)))))
			-> decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), !(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return  (!(_STD move(x) == _STD move(right.x))) || (!(_STD move(y) == _STD move(right.y))) || (!(_STD move(z) == _STD move(right.z)));
		}

		void operator=(ResultPacker_3<E>) = delete;

	};

	template<class E>
	struct ResultPacker_4 {
		E x,y,z,w;

		template<class T>
		EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
			EUCVECTORINLINE auto operator+(ResultPacker_4<T>&& pack) const noexcept(noexcept(ResultPacker_4<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)}))
			-> decltype(ResultPacker_4<meta::no_ref<decltype(_STD move(x) + _STD move(pack.x))>>{_STD move(x) + _STD move(pack.x)})  {
			return { _STD move(x) + _STD move(pack.x), _STD move(y) + _STD move(pack.y), _STD move(z) + _STD move(pack.z), _STD move(w) + _STD move(pack.w) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
			EUCVECTORINLINE auto operator-(ResultPacker_4<T>&& pack) const noexcept(noexcept(ResultPacker_4<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}))
			-> decltype(ResultPacker_4<meta::no_ref<decltype(_STD move(x) - _STD move(pack.x))>>{_STD move(x) - _STD move(pack.x)}) {
			return { _STD move(x) - _STD move(pack.x), _STD move(y) - _STD move(pack.y), _STD move(z) - _STD move(pack.z), _STD move(w) - _STD move(pack.w) };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE auto operator*(T&& scl) noexcept(noexcept(ResultPacker_4<meta::no_ref<decltype(x * _STD forward<T>(scl))>>{x * _STD forward<T>(scl)}))
			-> decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_4<meta::no_ref<decltype(x * _STD forward<T>(scl))>>{x * _STD forward<T>(scl)}) {
			return { x * scl,y * scl,z * scl,w * scl };
		}

		template<class T>
		friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
			EUCVECTORINLINE static auto operator*(T&& scl, ResultPacker_4<E>&& right) noexcept(noexcept(ResultPacker_4<meta::no_ref<decltype(right.x* scl)>>{right.x* scl}))
			-> decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_4<meta::no_ref<decltype(right.x* scl)>>{right.x* scl})  {
			return { right.x * scl,right.y * scl,right.z * scl,right.w * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
			EUCVECTORINLINE auto operator/(T&& scl) noexcept(noexcept(ResultPacker_4<meta::no_ref<decltype(x / _STD forward<T>(scl))>>{x / _STD forward<T>(scl)}))
			-> decltype(meta::cnd_type<!_STD is_base_of_v<meta::evd_euc_vec, meta::no_ref<T>>>(), ResultPacker_4<meta::no_ref<decltype(x / _STD forward<T>(scl))>>{x / _STD forward<T>(scl)}) {
			return { x * scl,y * scl,z * scl,w * scl };
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator==(ResultPacker_4<T>&& right) noexcept(noexcept(bool(_STD move(x) == _STD move(right.x))))
			-> decltype(bool(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) == _STD move(right.x)) && (_STD move(y) == _STD move(right.y)) && (_STD move(z) == _STD move(right.z)) && (_STD move(w) == _STD move(right.w));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_4<T>&& right) noexcept(noexcept(bool(_STD move(x) != _STD move(right.x))))
			-> decltype(bool(_STD move(x) != _STD move(right.x)), _STD declval<bool>()) {
			return (_STD move(x) != _STD move(right.x)) || (_STD move(y) != _STD move(right.y)) || (_STD move(z) != _STD move(right.z)) || (_STD move(w) != _STD move(right.w));
		}

		template<class T>
		EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
			EUCVECTORINLINE auto operator!=(ResultPacker_4<T>&& right) const noexcept(noexcept(bool(!(_STD move(x) == _STD move(right.x)))))
			-> decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), !(_STD move(x) == _STD move(right.x)), _STD declval<bool>()) {
			return  (!(_STD move(x) == _STD move(right.x))) || (!(_STD move(y) == _STD move(right.y))) || (!(_STD move(z) == _STD move(right.z))) || (!(_STD move(w) == _STD move(right.w)));
		}

		void operator=(ResultPacker_4<E>) = delete;

	};

};

/*
	Vector 1
*/
template<class E = float>
struct EuclideanVector1
	: protected meta::evd_euc_vec {
protected:

	static constexpr size_t EucD = 1;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using ElemType = E;
	using EucVector = EuclideanVector1<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_1<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanVector1;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

protected:

	ElemType x_;

public:

	/*
		Constructors.
	*/
	EuclideanVector1() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: x_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanVector1(const EuclideanVector1& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: x_(vector.x_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanVector1(EuclideanVector1&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: x_(_STD move(vector.x_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanVector1(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: x_(_STD move(pack.x))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	explicit EuclideanVector1(T&& val) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: x_(_STD forward<T>(val))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanVector1<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ + vector.x_)>>{x_ + vector.x_}))
		->decltype(Packer<meta::no_ref<decltype(x_ + vector.x_)>>{x_ + vector.x_}) {
		return { x_ + vector.x_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanVector1<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ + _STD move(vector.x_))>>{x_ + _STD move(vector.x_)}))
		->decltype(Packer<meta::no_ref<decltype(x_ + _STD move(vector.x_))>>{x_ + _STD move(vector.x_)}) {
		return { x_ + _STD move(vector.x_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ + _STD move(pack.x))>>{x_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(x_ + _STD move(pack.x))>>{x_ + _STD move(pack.x)}) {
		return { x_ + _STD move(pack.x) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.x_)>>{_STD move(pack.x) + vector.x_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.x_)>>{_STD move(pack.x) + vector.x_}) {
		return { _STD move(pack.x) + vector.x_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.x_))>>{_STD move(pack.x) + _STD move(vector.x_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.x_))>>{_STD move(pack.x) + _STD move(vector.x_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored.If you intend to modify the lvalue, use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanVector1<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ - vector.x_)>>{x_ - vector.x_}))
		->decltype(Packer<meta::no_ref<decltype(x_ - vector.x_)>>{x_ - vector.x_}) {
		return { x_ - vector.x_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored.If you intend to modify the lvalue, use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ - _STD move(pack.x))>>{x_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(x_ - _STD move(pack.x))>>{x_ - _STD move(pack.x)}) {
		return { x_ - _STD move(pack.x) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored.If you intend to modify the lvalue, use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanVector1<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ - _STD move(vector.x_))>>{x_ - _STD move(vector.x_)}))
		->decltype(Packer<meta::no_ref<decltype(x_ - _STD move(vector.x_))>>{x_ - _STD move(vector.x_)}) {
		return { x_ - _STD move(vector.x_) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored.If you intend to modify the lvalue, use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.x_)>>{_STD move(pack.x) - vector.x_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.x_)>>{_STD move(pack.x) - vector.x_}) {
		return { _STD move(pack.x) - vector.x_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored.If you intend to modify the lvalue, use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.x_))>>{_STD move(pack.x) - _STD move(vector.x_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.x_))>>{_STD move(pack.x) - _STD move(vector.x_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_* scl)>>{x_* scl}))
		->decltype(Packer<meta::no_ref<decltype(x_* scl)>>{x_* scl}) {
		return { x_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.x_* scl)>>{right.x_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.x_* scl)>>{right.x_* scl}) {
		return { right.x_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.x_)* scl)>>{_STD move(right.x_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.x_)* scl)>>{_STD move(right.x_)* scl}) {
		return { _STD move(right.x_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(x_ / scl)>>{x_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(x_ / scl)>>{x_ / scl}) {
		return { x_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanVector1<T>& right) const noexcept(noexcept(bool(x_ == right.x_)))
		-> meta::no_ref<decltype(bool(x_ == right.x_), _STD declval<bool>())> {
		return x_ == right.x_;
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanVector1<T>&& right) const noexcept(noexcept(bool(x_ == _STD move(right.x_))))
		-> meta::no_ref<decltype(bool(x_ == _STD move(right.x_)), _STD declval<bool>())> {
		return x_ == _STD move(right.x_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(x_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(x_ == _STD move(right.x)), _STD declval<bool>())> {
		return x_ == _STD move(right.x);
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.x_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.x_ == _STD move(left.x)), _STD declval<bool>())> {
		return right.x_ == _STD move(left.x);
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.x_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.x_) == _STD move(left.x)), _STD declval<bool>())> {
		return _STD move(right.x_) == _STD move(left.x);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanVector1<T>& right) const noexcept(noexcept(bool(x_ != right.x_)))
		-> meta::no_ref<decltype(bool(x_ != right.x_), _STD declval<bool>())> {
		return x_ != right.x_;
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanVector1<T>&& right) const noexcept(noexcept(bool(x_ != _STD move(right.x_))))
		-> meta::no_ref<decltype(bool(x_ != _STD move(right.x_)), _STD declval<bool>())> {
		return x_ != _STD move(right.x_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(x_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(x_ != _STD move(right.x)), _STD declval<bool>())> {
		return x_ != _STD move(right.x);
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.x_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.x_ != _STD move(left.x)), _STD declval<bool>())> {
		return right.x_ != _STD move(left.x);
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.x_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.x_) != _STD move(left.x)), _STD declval<bool>())> {
		return _STD move(right.x_) != _STD move(left.x);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanVector1<T>& right) const noexcept(noexcept(bool(!(x_ == right.x_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(x_ == right.x_)), _STD declval<bool>())> {
		return !(x_ == right.x_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanVector1<T>&& right) const noexcept(noexcept(bool(!(x_ == _STD move(right.x_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(x_ == _STD move(right.x_))), _STD declval<bool>())> {
		return !(x_ == _STD move(right.x_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(x_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(x_ == _STD move(right.x))), _STD declval<bool>())> {
		return !(x_ == _STD move(right.x));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.x_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.x_ == _STD move(left.x))), _STD declval<bool>())> {
		return !(right.x_ == _STD move(left.x));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.x_)) == _STD move(left.x)), _STD declval<bool>())> {
		return !(_STD move(right.x_) == _STD move(left.x));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanVector1<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanVector1<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanVector1<T>() * -1))
		->decltype(EuclideanVector1<T>() * -1) {
		return { x_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanVector1<T>& vector) & noexcept(noexcept(x_ = vector.x_))
		-> decltype(x_ = vector.x_, _STD declval<LRefEucVector>()) {
		x_ = vector.x_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanVector1<T>&& vector) & noexcept(noexcept(x_ = _STD move(vector.x_)))
		-> decltype(x_ = _STD move(vector.x_), _STD declval<LRefEucVector>()) {
		x_ = _STD move(vector.x_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(x_ = _STD move(pack.x)))
		-> decltype(x_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = _STD move(pack.x);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanVector1<T>& vector) & noexcept(noexcept(x_ += vector.x_))
		-> decltype(x_ += vector.x_, _STD declval<LRefEucVector>()) {
		x_ += vector.x_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanVector1<T>&& vector) & noexcept(noexcept(x_ += _STD move(vector.x_)))
		-> decltype(x_ += _STD move(vector.x_), _STD declval<LRefEucVector>()) {
		x_ += _STD move(vector.x_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(x_ += _STD move(pack.x)))
		-> decltype(x_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ += _STD move(pack.x);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanVector1<T>& vector) & noexcept(noexcept(x_ = x_ + vector.x_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, x_ = x_ + vector.x_, _STD declval<LRefEucVector>()) {
		x_ = x_ + vector.x_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanVector1<T>&& vector) & noexcept(noexcept(x_ = x_ + _STD move(vector.x_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, x_ = x_ + _STD move(vector.x_), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(vector.x_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(x_ = x_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, x_ = x_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(pack.x);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanVector1<T>& vector) & noexcept(noexcept(x_ -= vector.x_))
		-> decltype(x_ -= vector.x_, _STD declval<LRefEucVector>()) {
		x_ -= vector.x_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanVector1<T>&& vector) & noexcept(noexcept(x_ -= _STD move(vector.x_)))
		-> decltype(x_ -= _STD move(vector.x_), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(vector.x_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(x_ -= _STD move(pack.x)))
		-> decltype(x_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(pack.x);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanVector1<T>& vector) & noexcept(noexcept(x_ = x_ - vector.x_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, x_ = x_ - vector.x_, _STD declval<LRefEucVector>()) {
		x_ = x_ - vector.x_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanVector1<T>&& vector) & noexcept(noexcept(x_ = x_ - _STD move(vector.x_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, x_ = x_ - _STD move(vector.x_), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(vector.x_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(x_ = x_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, x_ = x_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(pack.x);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(x_ *= scl))
		-> decltype(x_ *= scl, _STD declval<LRefEucVector>()) {
		x_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(x_ = x_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, x_ = x_ * scl, _STD declval<LRefEucVector>()) {
		x_ = x_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(x_ /= scl))
		-> decltype(x_ /= scl, _STD declval<LRefEucVector>()) {
		x_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(x_ = x_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, x_ = x_ / scl, _STD declval<LRefEucVector>()) {
		x_ = x_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return x_; }
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		x_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class T>
	EUCVECTORINLINE auto set(T&& val) noexcept(noexcept(x_ = _STD forward<T>(val)))
		-> meta::is_type_t<decltype(x_ = _STD forward<T>(val)), void> {
		x_ = _STD forward<T>(val);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanVector1<T>& vector) const noexcept(noexcept(decltype(x_* vector.x_)(x_* vector.x_)))
		-> decltype(decltype(x_* vector.x_)(x_* vector.x_)) {
		return x_ * vector.x_;
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanVector1<T>&& vector) const noexcept(noexcept(decltype(x_* _STD move(vector.x_))(x_* _STD move(vector.x_))))
		-> decltype(decltype(x_* _STD move(vector.x_))(x_* _STD move(vector.x_))) {
		return x_ * _STD move(vector.x_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(decltype(x_* _STD move(pack.x))(x_* _STD move(pack.x))))
		-> decltype(decltype(x_* _STD move(pack.x))(x_* _STD move(pack.x))) {
		return x_ * _STD move(pack.x);
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return x_ * x_;
	}
	/*
		@brief

			Calculate the norm.

			out = ROOT((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return x_;
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{x_ / eucnorm<T>()}))
		->decltype(Packer<T>{x_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { x_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		auto&& norm = eucnorm<T>();
		x_ /= norm;
		return *this;
	}

};

/*
	Recursion Vector.
*/

/*
	D2.
*/
template<class E = float>
struct EuclideanRecVector2
	: protected EuclideanVector1<E> {
protected:

	static constexpr size_t EucD = 2;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using MX = EuclideanVector1<E>;
	using ElemType = E;
	using EucVector = EuclideanRecVector2<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_2<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanRecVector2;
	template<class FE>
	friend struct EuclideanCmplVector2;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

protected:

	ElemType y_;

public:

	/*
		Constructors.
	*/
	EuclideanRecVector2() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: MX()
		, y_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanRecVector2(const EuclideanRecVector2& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: MX(vector.x_)
		, y_(vector.y_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanRecVector2(EuclideanRecVector2&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: MX(_STD move(vector.x_))
		, y_(_STD move(vector.y_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanRecVector2(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: MX(_STD move(pack.x))
		, y_(_STD move(pack.y))
	{}

	template<class X, class Y, meta::if_t<meta::is_constructible_anynum_param_v<ElemType,X,Y>> = 0>
	EuclideanRecVector2(X&& x, Y&& y) noexcept(noexcept(MX(_STD forward<X>(x))) && _STD is_nothrow_constructible_v<ElemType, Y>)
		: MX(_STD forward<X>(x))
		, y_(_STD forward<Y>(y))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanRecVector2<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + vector.y_)>>{y_ + vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(y_ + vector.y_)>>{y_ + vector.y_}) {
		return { MX::x_ + vector.x_, y_ + vector.y_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanRecVector2<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + _STD move(vector.y_))>>{y_ + _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(y_ + _STD move(vector.y_))>>{y_ + _STD move(vector.y_)}) {
		return { MX::x_ + _STD move(vector.x_), y_ + _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + _STD move(pack.x))>>{y_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(y_ + _STD move(pack.x))>>{y_ + _STD move(pack.x)}) {
		return { MX::x_ + _STD move(pack.x), y_ + _STD move(pack.y) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.y_)>>{_STD move(pack.x) + vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.y_)>>{_STD move(pack.x) + vector.y_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.y_))>>{_STD move(pack.x) + _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.y_))>>{_STD move(pack.x) + _STD move(vector.y_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_),_STD move(pack.y) + _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanRecVector2<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - vector.y_)>>{y_ - vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(y_ - vector.y_)>>{y_ - vector.y_}) {
		return { MX::x_ - vector.x_, y_ - vector.y_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanRecVector2<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - _STD move(vector.y_))>>{y_ - _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(y_ - _STD move(vector.y_))>>{y_ - _STD move(vector.y_)}) {
		return { MX::x_ - _STD move(vector.x_), y_ - _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - _STD move(pack.x))>>{y_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(y_ - _STD move(pack.x))>>{y_ - _STD move(pack.x)}) {
		return { MX::x_ - _STD move(pack.x), y_ - _STD move(pack.y) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.y_)>>{_STD move(pack.x) - vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.y_)>>{_STD move(pack.x) - vector.y_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.y_))>>{_STD move(pack.x) - _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.y_))>>{_STD move(pack.x) - _STD move(vector.y_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_),_STD move(pack.y) - _STD move(vector.y_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_* scl)>>{y_* scl}))
		->decltype(Packer<meta::no_ref<decltype(y_* scl)>>{y_* scl}) {
		return { MX::x_ * scl, y_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.y_* scl)>>{right.y_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.y_* scl)>>{right.y_* scl}) {
		return { right.x_ * scl, right.y_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.y_)* scl)>>{_STD move(right.y_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.y_)* scl)>>{_STD move(right.y_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ / scl)>>{y_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(y_ / scl)>>{y_ / scl}) {
		return { MX::x_ / scl, y_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanRecVector2<T>& right) const noexcept(noexcept(bool(y_ == right.y_)))
		-> meta::no_ref<decltype(bool(y_ == right.y_))> {
		return (MX::x_ == right.x_) && (y_ == right.y_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanRecVector2<T>&& right) const noexcept(noexcept(bool(y_ == _STD move(right.y_))))
		-> meta::no_ref<decltype(bool(y_ == _STD move(right.y_)))> {
		return (MX::x_ == _STD move(right.x_)) && (y_ == _STD move(right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(y_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(y_ == _STD move(right.x)))> {
		return (MX::x_ == _STD move(right.x)) && (y_ == _STD move(right.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.y_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.y_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.y_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.y_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector2<T>& right) const noexcept(noexcept(bool(y_ != right.y_)))
		-> meta::no_ref<decltype(bool(y_ != right.y_))> {
		return (MX::x_ != right.x_) || (y_ != right.y_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector2<T>&& right) const noexcept(noexcept(bool(y_ != _STD move(right.y_))))
		-> meta::no_ref<decltype(bool(y_ != _STD move(right.y_)))> {
		return (MX::x_ != _STD move(right.x_)) || (y_ != _STD move(right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(y_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(y_ != _STD move(right.x)), _STD declval<bool>())> {
		return (MX::x_ != _STD move(right.x)) || (y_ != _STD move(right.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.y_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.y_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.y_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.y_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector2<T>& right) const noexcept(noexcept(bool(!(y_ == right.y_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(y_ == right.y_)), _STD declval<bool>())> {
		return (!(MX::x_ == right.x_)) || (!(y_ == right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector2<T>&& right) const noexcept(noexcept(bool(!(y_ == _STD move(right.y_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(y_ == _STD move(right.y_))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x_))) || (!(y_ == _STD move(right.y_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(y_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(y_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x))) || (!(y_ == _STD move(right.y)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.y_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.y_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.y_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanRecVector2<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanRecVector2<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanRecVector2<T>() * -1))
		->decltype(EuclideanRecVector2<T>() * -1) {
		return { MX::x_ * -1, y_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanRecVector2<T>& vector) & noexcept(noexcept(y_ = vector.y_))
		-> decltype(y_ = vector.y_, _STD declval<LRefEucVector>()) {
		MX::x_ = vector.x_;
		y_ = vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanRecVector2<T>&& vector) & noexcept(noexcept(y_ = _STD move(vector.y_)))
		-> decltype(y_ = _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(vector.x_);
		y_ = _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(y_ = _STD move(pack.x)))
		-> decltype(y_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(pack.x);
		y_ = _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector2<T>& vector) & noexcept(noexcept(y_ += vector.y_))
		-> decltype(y_ += vector.y_, _STD declval<LRefEucVector>()) {
		MX::x_ += vector.x_;
		y_ += vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector2<T>&& vector) & noexcept(noexcept(y_ += _STD move(vector.y_)))
		-> decltype(y_ += _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(vector.x_);
		y_ += _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(y_ += _STD move(pack.x)))
		-> decltype(y_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(pack.x);
		y_ += _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector2<T>& vector) & noexcept(noexcept(y_ = y_ + vector.y_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, y_ = y_ + vector.y_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + vector.x_;
		y_ = y_ + vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector2<T>&& vector) & noexcept(noexcept(y_ = y_ + _STD move(vector.y_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, y_ = y_ + _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(vector.x_);
		y_ = y_ + _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(y_ = y_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, y_ = y_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(pack.x);
		y_ = y_ + _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector2<T>& vector) & noexcept(noexcept(y_ -= vector.y_))
		-> decltype(y_ -= vector.y_, _STD declval<LRefEucVector>()) {
		MX::x_ -= vector.x_;
		y_ -= vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector2<T>&& vector) & noexcept(noexcept(y_ -= _STD move(vector.y_)))
		-> decltype(y_ -= _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(vector.x_);
		y_ -= _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(y_ -= _STD move(pack.x)))
		-> decltype(y_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(pack.x);
		y_ -= _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector2<T>& vector) & noexcept(noexcept(y_ = y_ - vector.y_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, y_ = y_ - vector.y_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - vector.x_;
		y_ = y_ - vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector2<T>&& vector) & noexcept(noexcept(y_ = y_ - _STD move(vector.y_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, y_ = y_ - _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(vector.x_);
		y_ = y_ - _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(y_ = y_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, y_ = y_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(pack.x);
		y_ = y_ - _STD move(pack.y);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(y_ *= scl))
		-> decltype(y_ *= scl, _STD declval<LRefEucVector>()) {
		MX::x_ *= scl;
		y_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(y_ = y_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, y_ = y_ * scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ * scl;
		y_ = y_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(y_ /= scl))
		-> decltype(y_ /= scl, _STD declval<LRefEucVector>()) {
		MX::x_ /= scl;
		y_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(y_ = y_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, y_ = y_ / scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ / scl;
		y_ = y_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xx() const noexcept { return { MX::x_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yx() const noexcept { return { y_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yy() const noexcept { return { y_,y_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { MX::x_,y_ }; }
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		MX::x_ = 0;
		y_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y>
	EUCVECTORINLINE auto set(X&& x, Y&& y) noexcept(noexcept(MX::set(_STD forward<X>(x))) && noexcept(y_ = _STD forward<Y>(y)))
		-> meta::is_type_t<decltype(MX::set(_STD forward<X>(x)), y_ = _STD forward<Y>(y)), void> {
		MX::x_ = _STD forward<X>(x);
		y_ = _STD forward<Y>(y);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanRecVector2<T>& vector) const noexcept(noexcept(MX::dot(vector) + decltype(y_ * vector.y_)(y_ * vector.y_)))
		-> decltype(MX::dot(vector) + decltype(y_ * vector.y_)(y_ * vector.y_)) {
		return (MX::x_ * vector.x_) + (y_ * vector.y_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanRecVector2<T>&& vector) const noexcept(noexcept(MX::dot(_STD move(vector)) + decltype(y_ * _STD move(vector.y_))(y_ * _STD move(vector.y_))))
		-> decltype(MX::dot(_STD move(vector)) + decltype(y_ * _STD move(vector.y_))(y_ * _STD move(vector.y_))) {
		return (MX::x_ * _STD move(vector.x_)) + (y_ * _STD move(vector.y_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(MX::dot(_STD move(pack)) + decltype(y_* _STD move(pack.x))(y_* _STD move(pack.x))))
		-> decltype(MX::dot(_STD move(pack)) + decltype(y_ * _STD move(pack.x))(y_ * _STD move(pack.x))) {
		return (MX::x_ * _STD move(pack.x)) + (y_ * _STD move(pack.y));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return (MX::x_* MX::x_) + (y_ * y_);
	}
	/*
		@brief

			Calculate the norm.

			out = ROOT((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt((MX::x_ * MX::x_) + (y_ * y_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{y_ / eucnorm<T>()}))
		->decltype(Packer<T>{y_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { MX::x_ / norm, y_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}

};

/*
	D3.
*/
template<class E = float>
struct EuclideanRecVector3
	: protected EuclideanRecVector2<E> {
protected:

	static constexpr size_t EucD = 3;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using MX = EuclideanVector1<E>;
	using MXY = EuclideanRecVector2<E>;
	using ElemType = E;
	using EucVector = EuclideanRecVector3<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_3<R>;
	template<class R>
	using SubPacker2 = detail::ResultPacker_2<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanRecVector3;
	template<class FE>
	friend struct EuclideanCmplVector3;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

protected:

	ElemType z_;

public:

	/*
		Constructors.
	*/
	EuclideanRecVector3() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: MXY()
		, z_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanRecVector3(const EuclideanRecVector3& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: MXY(vector)
		, z_(vector.z_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanRecVector3(EuclideanRecVector3&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: MXY(_STD move(vector))
		, z_(_STD move(vector.z_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanRecVector3(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: MXY(_STD move(pack.x), _STD move(pack.y))
		, z_(_STD move(pack.z))
	{}

	template<class X, class Y,class Z, meta::if_t<meta::is_constructible_anynum_param_v<ElemType, X, Y, Z>> = 0>
	EuclideanRecVector3(X&& x, Y&& y, Z&& z) noexcept(noexcept(MXY(_STD forward<X>(x), _STD forward<Y>(y))) && _STD is_nothrow_constructible_v<ElemType, Z>)
		: MXY(_STD forward<X>(x), _STD forward<Y>(y))
		, z_(_STD forward<Z>(z))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanRecVector3<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + vector.z_)>>{z_ + vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(z_ + vector.z_)>>{z_ + vector.z_}) {
		return { MX::x_ + vector.x_, MXY::y_ + vector.y_, z_ + vector.z_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanRecVector3<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + _STD move(vector.z_))>>{z_ + _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(z_ + _STD move(vector.z_))>>{z_ + _STD move(vector.z_)}) {
		return { MX::x_ + _STD move(vector.x_), MXY::y_ + _STD move(vector.y_), z_ + _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + _STD move(pack.x))>>{z_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(z_ + _STD move(pack.x))>>{z_ + _STD move(pack.x)}) {
		return { MX::x_ + _STD move(pack.x), MXY::y_ + _STD move(pack.y), z_ + _STD move(pack.z) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.z_)>>{_STD move(pack.x) + vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.z_)>>{_STD move(pack.x) + vector.z_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ , _STD move(pack.z) + vector.z_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.z_))>>{_STD move(pack.x) + _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.z_))>>{_STD move(pack.x) + _STD move(vector.z_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_), _STD move(pack.y) + _STD move(vector.y_), _STD move(pack.z) + _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanRecVector3<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - vector.z_)>>{z_ - vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(z_ - vector.z_)>>{z_ - vector.z_}) {
		return { MX::x_ - vector.x_, MXY::y_ - vector.y_, z_ - vector.z_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanRecVector3<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - _STD move(vector.z_))>>{z_ - _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(z_ - _STD move(vector.z_))>>{z_ - _STD move(vector.z_)}) {
		return { MX::x_ - _STD move(vector.x_), MXY::y_ - _STD move(vector.y_), z_ - _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - _STD move(pack.x))>>{z_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(z_ - _STD move(pack.x))>>{z_ - _STD move(pack.x)}) {
		return { MX::x_ - _STD move(pack.x), MXY::y_ - _STD move(pack.y), z_ - _STD move(pack.z) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.z_)>>{_STD move(pack.x) - vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.z_)>>{_STD move(pack.x) - vector.z_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ , _STD move(pack.z) - vector.z_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.z_))>>{_STD move(pack.x) - _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.z_))>>{_STD move(pack.x) - _STD move(vector.z_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_), _STD move(pack.y) - _STD move(vector.y_), _STD move(pack.z) - _STD move(vector.z_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_* scl)>>{z_* scl}))
		->decltype(Packer<meta::no_ref<decltype(z_* scl)>>{z_* scl}) {
		return { MX::x_ * scl, MXY::y_ * scl, z_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.z_* scl)>>{right.z_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.z_* scl)>>{right.z_* scl}) {
		return { right.x_ * scl,  right.y_ * scl, right.z_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.z_)* scl)>>{_STD move(right.z_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.z_)* scl)>>{_STD move(right.z_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl, _STD move(right.z_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ / scl)>>{z_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(z_ / scl)>>{z_ / scl}) {
		return { MX::x_ / scl, MXY::y_ / scl, z_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanRecVector3<T>& right) const noexcept(noexcept(bool(z_ == right.z_)))
		-> meta::no_ref<decltype(bool(z_ == right.z_))> {
		return (MX::x_ == right.x_) && (MXY::y_ == right.y_) && (z_ == right.z_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanRecVector3<T>&& right) const noexcept(noexcept(bool(z_ == _STD move(right.z_))))
		-> meta::no_ref<decltype(bool(z_ == _STD move(right.z_)))> {
		return (MX::x_ == _STD move(right.x_)) && (MXY::y_ == _STD move(right.y_)) && (z_ == _STD move(right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(z_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(z_ == _STD move(right.x)))> {
		return (MX::x_ == _STD move(right.x)) && (MXY::y_ == _STD move(right.y)) && (z_ == _STD move(right.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.z_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.z_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y)) && (right.z_ == _STD move(left.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.z_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.z_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y)) && (_STD move(right.z_) == _STD move(left.z));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector3<T>& right) const noexcept(noexcept(bool(z_ != right.z_)))
		-> meta::no_ref<decltype(bool(z_ != right.z_))> {
		return (MX::x_ != right.x_) || (MXY::y_ != right.y_) || (z_ != right.z_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector3<T>&& right) const noexcept(noexcept(bool(z_ != _STD move(right.z_))))
		-> meta::no_ref<decltype(bool(z_ != _STD move(right.z_)))> {
		return (MX::x_ != _STD move(right.x_)) || (MXY::y_ != _STD move(right.y_)) || (z_ != _STD move(right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(z_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(z_ != _STD move(right.x)), _STD declval<bool>())> {
		return (MX::x_ != _STD move(right.x)) || (MXY::y_ != _STD move(right.y)) || (z_ != _STD move(right.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.z_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.z_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y)) || (right.z_ != _STD move(left.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.z_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.z_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y)) || (_STD move(right.z_) != _STD move(left.z));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector3<T>& right) const noexcept(noexcept(bool(!(z_ == right.z_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(z_ == right.z_)), _STD declval<bool>())> {
		return (!(MX::x_ == right.x_)) || (!(MXY::y_ == right.y_)) || (!(z_ == right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector3<T>&& right) const noexcept(noexcept(bool(!(z_ == _STD move(right.z_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(z_ == _STD move(right.z_))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x_))) || (!(MXY::y_ == _STD move(right.y_))) || (!(z_ == _STD move(right.z_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(z_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(z_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x))) || (!(MXY::y_ == _STD move(right.y))) || (!(z_ == _STD move(right.z)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.z_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.z_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y))) || (!(right.z_ == _STD move(left.z)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.z_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y))) || (!(_STD move(right.z_) == _STD move(left.z)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanRecVector3<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanRecVector3<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanRecVector3<T>() * -1))
		->decltype(EuclideanRecVector3<T>() * -1) {
		return { MX::x_ * -1, MXY::y_ * -1, z_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanRecVector3<T>& vector) & noexcept(noexcept(z_ = vector.z_))
		-> decltype(z_ = vector.z_, _STD declval<LRefEucVector>()) {
		MX::x_ = vector.x_;
		MXY::y_ = vector.y_;
		z_ = vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanRecVector3<T>&& vector) & noexcept(noexcept(z_ = _STD move(vector.z_)))
		-> decltype(z_ = _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(vector.x_);
		MXY::y_ = _STD move(vector.y_);
		z_ = _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(z_ = _STD move(pack.x)))
		-> decltype(z_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(pack.x);
		MXY::y_ = _STD move(pack.y);
		z_ = _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector3<T>& vector) & noexcept(noexcept(z_ += vector.z_))
		-> decltype(z_ += vector.z_, _STD declval<LRefEucVector>()) {
		MX::x_ += vector.x_;
		MXY::y_ += vector.y_;
		z_ += vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector3<T>&& vector) & noexcept(noexcept(z_ += _STD move(vector.z_)))
		-> decltype(z_ += _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(vector.x_);
		MXY::y_ += _STD move(vector.y_);
		z_ += _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(z_ += _STD move(pack.x)))
		-> decltype(z_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(pack.x);
		MXY::y_ += _STD move(pack.y);
		z_ += _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector3<T>& vector) & noexcept(noexcept(z_ = z_ + vector.z_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, z_ = z_ + vector.z_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + vector.x_;
		MXY::y_ = MXY::y_ + vector.y_;
		z_ = z_ + vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector3<T>&& vector) & noexcept(noexcept(z_ = z_ + _STD move(vector.z_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, z_ = z_ + _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(vector.x_);
		MXY::y_ = MXY::y_ + _STD move(vector.y_);
		z_ = z_ + _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(z_ = z_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, z_ = z_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(pack.x);
		MXY::y_ = MXY::y_ + _STD move(pack.y);
		z_ = z_ + _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector3<T>& vector) & noexcept(noexcept(z_ -= vector.z_))
		-> decltype(z_ -= vector.z_, _STD declval<LRefEucVector>()) {
		MX::x_ -= vector.x_;
		MXY::y_ -= vector.y_;
		z_ -= vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector3<T>&& vector) & noexcept(noexcept(z_ -= _STD move(vector.z_)))
		-> decltype(z_ -= _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(vector.x_);
		MXY::y_ -= _STD move(vector.y_);
		z_ -= _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(z_ -= _STD move(pack.x)))
		-> decltype(z_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(pack.x);
		MXY::y_ -= _STD move(pack.y);
		z_ -= _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector3<T>& vector) & noexcept(noexcept(z_ = z_ - vector.z_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, z_ = z_ - vector.z_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - vector.x_;
		MXY::y_ = MXY::y_ - vector.y_;
		z_ = z_ - vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector3<T>&& vector) & noexcept(noexcept(z_ = z_ - _STD move(vector.z_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, z_ = z_ - _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(vector.x_);
		MXY::y_ = MXY::y_ - _STD move(vector.y_);
		z_ = z_ - _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(z_ = z_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, z_ = z_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(pack.x);
		MXY::y_ = MXY::y_ - _STD move(pack.y);
		z_ = z_ - _STD move(pack.z);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(z_ *= scl))
		-> decltype(z_ *= scl, _STD declval<LRefEucVector>()) {
		MX::x_ *= scl;
		MXY::y_ *= scl;
		z_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(z_ = z_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, z_ = z_ * scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ * scl;
		MXY::y_ = MXY::y_ * scl;
		z_ = z_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(z_ /= scl))
		-> decltype(z_ /= scl, _STD declval<LRefEucVector>()) {
		MX::x_ /= scl;
		MXY::y_ /= scl;
		z_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(z_ = z_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, z_ = z_ / scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ / scl;
		MXY::y_ = MXY::y_ / scl;
		z_ = z_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return MXY::y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return  MXY::y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType z() noexcept { return z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType z() const noexcept { return z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr MXY& xy() noexcept { return *this; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr const MXY& xy() const noexcept { return *this; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xx() const noexcept { return { MX::x_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xz() const noexcept { return { MX::x_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yx() const noexcept { return { MXY::y_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yy() const noexcept { return { MXY::y_,MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yz() const noexcept { return { MXY::y_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zx() const noexcept { return { z_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zy() const noexcept { return { z_,MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zz() const noexcept { return { z_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxx() const noexcept { return { MX::x_,MX::x_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxy() const noexcept { return { MX::x_,MX::x_,MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxz() const noexcept { return { MX::x_, MX::x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyx() const noexcept { return { MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyy() const noexcept { return { MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyz() const noexcept { return { MX::x_, MXY::y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzx() const noexcept { return { MX::x_, z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzy() const noexcept { return { MX::x_, z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzz() const noexcept { return { MX::x_, z_, z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxx() const noexcept { return { MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxy() const noexcept { return { MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxz() const noexcept { return { MXY::y_, MX::x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyx() const noexcept { return { MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyy() const noexcept { return { MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyz() const noexcept { return { MXY::y_, MXY::y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzx() const noexcept { return { MXY::y_, z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzy() const noexcept { return { MXY::y_, z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzz() const noexcept { return { MXY::y_, z_, z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxx() const noexcept { return { z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxy() const noexcept { return { z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxz() const noexcept { return { z_, MX::x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyx() const noexcept { return { z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyy() const noexcept { return { z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyz() const noexcept { return { z_, MXY::y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzx() const noexcept { return { z_, z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzy() const noexcept { return { z_, z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzz() const noexcept { return { z_, z_, z_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { MX::x_,MXY::y_, z_ }; }
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		MX::x_ = 0;
		MXY::y_ = 0;
		z_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y, class Z>
	EUCVECTORINLINE auto set(X&& x, Y&& y, Z&& z) noexcept(noexcept(MXY::set(_STD forward<X>(x), _STD forward<Y>(y))) && noexcept(z_ = _STD forward<Z>(z)))
		-> meta::is_type_t<decltype(MXY::set(_STD forward<X>(x), _STD forward<Y>(y)), z_ = _STD forward<Z>(z)), void> {
		MX::x_ = _STD forward<X>(x);
		MXY::y_ = _STD forward<Y>(y);
		z_ = _STD forward<Z>(z);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanRecVector3<T>& vector) const noexcept(noexcept(MXY::dot(vector) + decltype(z_ * vector.z_)(z_ * vector.z_)))
		-> decltype(MXY::dot(vector) + decltype(z_ * vector.z_)(z_ * vector.z_)) {
		return (MX::x_ * vector.x_) + (MXY::y_ * vector.y_) + (z_ * vector.z_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanRecVector3<T>&& vector) const noexcept(noexcept(MXY::dot(_STD move(vector)) + decltype(z_ * _STD move(vector.z_))(z_ * _STD move(vector.z_))))
		-> decltype(MXY::dot(_STD move(vector)) + decltype(z_ * _STD move(vector.z_))(z_ * _STD move(vector.z_))) {
		return (MX::x_ * _STD move(vector.x_)) + (MXY::y_ * _STD move(vector.y_)) + (z_ * _STD move(vector.z_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(MXY::dot(_STD move(pack)) + decltype(z_ * _STD move(pack.x))(z_ * _STD move(pack.x))))
		-> decltype(MXY::dot(_STD move(pack)) + decltype(z_ * _STD move(pack.x))(z_ * _STD move(pack.x))) {
		return (MX::x_ * _STD move(pack.x)) + (MXY::y_ * _STD move(pack.y)) + (z_ * _STD move(pack.z));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return (MX::x_ * MX::x_) + (MXY::y_ * MXY::y_) + (z_ * z_);
	}
	/*
		@brief

			Calculate the norm.

			out = root((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt((MX::x_ * MX::x_) + (MXY::y_ * MXY::y_) + (z_ * z_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{z_ / eucnorm<T>()}))
		->decltype(Packer<T>{z_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { MX::x_ / norm, MXY::y_ / norm, z_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}
	/*
		@brief

			Calculate the cross product.
	
	*/
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(const EuclideanRecVector3<T>& right) const noexcept(
			noexcept(Packer<decltype(z_* right.z_ - z_ * right.z_)>{ z_* right.z_ - z_ * right.z_ }))
		-> decltype(Packer<decltype(z_* right.z_ - z_ * right.z_)>{ z_* right.z_ - z_ * right.z_ }) {
		return { MXY::y_ * right.z_ - z_ * right.y_, z_ * right.x_ - MX::x_ * right.z_, MX::x_ * right.y_ - MXY::y_ * right.x_ };
	}
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(EuclideanRecVector3<T>&& right) const noexcept(
			noexcept(Packer<decltype(z_* _STD move(right.z_) - z_ * _STD move(right.z_))>{ z_* _STD move(right.z_) - z_ * _STD move(right.z_) }))
		-> decltype(Packer<decltype(z_* _STD move(right.z_) - z_ * _STD move(right.z_))>{ z_* _STD move(right.z_) - z_ * _STD move(right.z_) }) {
		return { MXY::y_ * _STD move(right.z_) - z_ * _STD move(right.y_), z_ * _STD move(right.x_) - MX::x_ * _STD move(right.z_), MX::x_ * _STD move(right.y_) - MXY::y_ * _STD move(right.x_) };
	}
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(RRefPacker<T> right) const noexcept(
			noexcept(Packer<decltype(z_* _STD move(right.z) - z_ * _STD move(right.z))>{ z_* _STD move(right.z) - z_ * _STD move(right.z) }))
		-> decltype(Packer<decltype(z_* _STD move(right.z) - z_ * _STD move(right.z))>{ z_* _STD move(right.z) - z_ * _STD move(right.z) }) {
		return { MXY::y_ * _STD move(right.z) - z_ * _STD move(right.y), z_ * _STD move(right.x) - MX::x_ * _STD move(right.z), MX::x_ * _STD move(right.y) - MXY::y_ * _STD move(right.x) };
	}
};

/*
	D4.
*/
template<class E = float>
struct EuclideanRecVector4
	: protected EuclideanRecVector3<E> {
protected:

	static constexpr size_t EucD = 4;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using MX = EuclideanVector1<E>;
	using MXY = EuclideanRecVector2<E>;
	using MXYZ = EuclideanRecVector3<E>;
	using ElemType = E;
	using EucVector = EuclideanRecVector4<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_4<R>;
	template<class R>
	using SubPacker2 = detail::ResultPacker_2<R>;
	template<class R>
	using SubPacker3 = detail::ResultPacker_3<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanRecVector4;
	template<class FE>
	friend struct EuclideanCmplVector4;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

protected:

	ElemType w_;

public:

	/*
		Constructors.
	*/
	EuclideanRecVector4() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: MXYZ()
		, w_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanRecVector4(const EuclideanRecVector4& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: MXYZ(vector)
		, w_(vector.w_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanRecVector4(EuclideanRecVector4&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: MXYZ(_STD move(vector))
		, w_(_STD move(vector.w_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanRecVector4(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: MXYZ(_STD move(pack.x), _STD move(pack.y), _STD move(pack.z))
		, w_(_STD move(pack.w))
	{}

	template<class X, class Y, class Z, class W, meta::if_t<meta::is_constructible_anynum_param_v<ElemType, X, Y, Z, W>> = 0>
	EuclideanRecVector4(X&& x, Y&& y, Z&& z, W&& w) noexcept(noexcept(MXYZ(_STD forward<X>(x), _STD forward<Y>(y), _STD forward<Z>(z))) && _STD is_nothrow_constructible_v<ElemType, W>)
		: MXYZ(_STD forward<X>(x), _STD forward<Y>(y), _STD forward<Z>(z))
		, w_(_STD forward<W>(w))
	{}
	
	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanRecVector4<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + vector.w_)>>{w_ + vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(w_ + vector.w_)>>{w_ + vector.w_}) {
		return { MX::x_ + vector.x_, MXY::y_ + vector.y_, MXYZ::z_ + vector.z_, w_ + vector.w_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanRecVector4<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + _STD move(vector.w_))>>{w_ + _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(w_ + _STD move(vector.w_))>>{w_ + _STD move(vector.w_)}) {
		return { MX::x_ + _STD move(vector.x_), MXY::y_ + _STD move(vector.y_), MXYZ::z_ + _STD move(vector.z_), w_ + _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + _STD move(pack.x))>>{w_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(w_ + _STD move(pack.x))>>{w_ + _STD move(pack.x)}) {
		return { MX::x_ + _STD move(pack.x), MXY::y_ + _STD move(pack.y), MXYZ::z_ + _STD move(pack.z), w_ + _STD move(pack.w) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.w_)>>{_STD move(pack.x) + vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.w_)>>{_STD move(pack.x) + vector.w_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ , _STD move(pack.z) + vector.z_ , _STD move(pack.w) + vector.w_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.w_))>>{_STD move(pack.x) + _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.w_))>>{_STD move(pack.x) + _STD move(vector.w_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_), _STD move(pack.y) + _STD move(vector.y_),  _STD move(pack.z) + _STD move(vector.z_), _STD move(pack.w) + _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanRecVector4<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - vector.w_)>>{w_ - vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(w_ - vector.w_)>>{w_ - vector.w_}) {
		return { MX::x_ - vector.x_, MXY::y_ - vector.y_, MXYZ::z_ - vector.z_, w_ - vector.w_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanRecVector4<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - _STD move(vector.w_))>>{w_ - _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(w_ - _STD move(vector.w_))>>{w_ - _STD move(vector.w_)}) {
		return { MX::x_ - _STD move(vector.x_), MXY::y_ - _STD move(vector.y_), MXYZ::z_ - _STD move(vector.z_), w_ - _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - _STD move(pack.x))>>{w_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(w_ - _STD move(pack.x))>>{w_ - _STD move(pack.x)}) {
		return { MX::x_ - _STD move(pack.x), MXY::y_ - _STD move(pack.y), MXYZ::z_ - _STD move(pack.z), w_ - _STD move(pack.w) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.w_)>>{_STD move(pack.x) - vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.w_)>>{_STD move(pack.x) - vector.w_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ , _STD move(pack.z) - vector.z_ , _STD move(pack.w) - vector.w_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.w_))>>{_STD move(pack.x) - _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.w_))>>{_STD move(pack.x) - _STD move(vector.w_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_), _STD move(pack.y) - _STD move(vector.y_),  _STD move(pack.z) - _STD move(vector.z_), _STD move(pack.w) - _STD move(vector.w_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_* scl)>>{w_* scl}))
		->decltype(Packer<meta::no_ref<decltype(w_* scl)>>{w_* scl}) {
		return { MX::x_ * scl, MXY::y_ * scl, MXYZ::z_ * scl, w_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.w_* scl)>>{right.w_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.w_* scl)>>{right.w_* scl}) {
		return { right.x_ * scl,  right.y_ * scl, right.z_ * scl, right.w_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.w_)* scl)>>{_STD move(right.w_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.w_)* scl)>>{_STD move(right.w_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl, _STD move(right.z_) * scl, _STD move(right.w_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ / scl)>>{w_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(w_ / scl)>>{w_ / scl}) {
		return { MX::x_ / scl, MXY::y_ / scl,  MXYZ::z_ / scl, w_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanRecVector4<T>& right) const noexcept(noexcept(bool(w_ == right.w_)))
		-> meta::no_ref<decltype(bool(w_ == right.w_))> {
		return (MX::x_ == right.x_) && (MXY::y_ == right.y_) && (MXYZ::z_ == right.z_) && (w_ == right.w_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanRecVector4<T>&& right) const noexcept(noexcept(bool(w_ == _STD move(right.w_))))
		-> meta::no_ref<decltype(bool(w_ == _STD move(right.w_)))> {
		return (MX::x_ == _STD move(right.x_)) && (MXY::y_ == _STD move(right.y_)) && (MXYZ::z_ == _STD move(right.z_)) && (w_ == _STD move(right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(w_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(w_ == _STD move(right.x)))> {
		return (MX::x_ == _STD move(right.x)) && (MXY::y_ == _STD move(right.y)) && (MXYZ::z_ == _STD move(right.z)) && (w_ == _STD move(right.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.w_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.w_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y)) && (right.z_ == _STD move(left.z)) && (right.w_ == _STD move(left.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.w_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.w_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y)) && (_STD move(right.z_) == _STD move(left.z)) && (_STD move(right.w_) == _STD move(left.w));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector4<T>& right) const noexcept(noexcept(bool(w_ != right.w_)))
		-> meta::no_ref<decltype(bool(w_ != right.w_))> {
		return (MX::x_ != right.x_) || (MXY::y_ != right.y_) || (MXYZ::z_ != right.z_) || (w_ != right.w_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector4<T>&& right) const noexcept(noexcept(bool(w_ != _STD move(right.w_))))
		-> meta::no_ref<decltype(bool(w_ != _STD move(right.w_)))> {
		return (MX::x_ != _STD move(right.x_)) || (MXY::y_ != _STD move(right.y_)) || (MXYZ::z_ != _STD move(right.z_)) || (w_ != _STD move(right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(w_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(w_ != _STD move(right.x)), _STD declval<bool>())> {
		return (MX::x_ != _STD move(right.x)) || (MXY::y_ != _STD move(right.y)) || (MXYZ::z_ != _STD move(right.z)) || (w_ != _STD move(right.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.w_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.w_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y)) || (right.z_ != _STD move(left.z)) || (right.w_ != _STD move(left.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.w_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.w_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y)) || (_STD move(right.z_) != _STD move(left.z)) || (_STD move(right.w_) != _STD move(left.w));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanRecVector4<T>& right) const noexcept(noexcept(bool(!(w_ == right.w_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(w_ == right.w_)), _STD declval<bool>())> {
		return (!(MX::x_ == right.x_)) || (!(MXY::y_ == right.y_)) || (!(MXYZ::z_ == right.z_)) || (!(w_ == right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanRecVector4<T>&& right) const noexcept(noexcept(bool(!(w_ == _STD move(right.w_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(w_ == _STD move(right.w_))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x_))) || (!(MXY::y_ == _STD move(right.y_))) || (!(MXYZ::z_ == _STD move(right.z_))) || (!(w_ == _STD move(right.w_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(w_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(w_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(MX::x_ == _STD move(right.x))) || (!(MXY::y_ == _STD move(right.y))) || (!(MXYZ::z_ == _STD move(right.z))) || (!(w_ == _STD move(right.w)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.w_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.w_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y))) || (!(right.z_ == _STD move(left.z))) || (!(right.w_ == _STD move(left.w)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.w_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y))) || (!(_STD move(right.z_) == _STD move(left.z))) || (!(_STD move(right.w_) == _STD move(left.w)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanRecVector4<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanRecVector4<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanRecVector4<T>() * -1))
		->decltype(EuclideanRecVector4<T>() * -1) {
		return { MX::x_ * -1, MXY::y_ * -1,  MXYZ::z_ * -1, w_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanRecVector4<T>& vector) & noexcept(noexcept(w_ = vector.w_))
		-> decltype(w_ = vector.w_, _STD declval<LRefEucVector>()) {
		MX::x_ = vector.x_;
		MXY::y_ = vector.y_;
		MXYZ::z_ = vector.z_;
		w_ = vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanRecVector4<T>&& vector) & noexcept(noexcept(w_ = _STD move(vector.w_)))
		-> decltype(w_ = _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(vector.x_);
		MXY::y_ = _STD move(vector.y_);
		MXYZ::z_ = _STD move(vector.z_);
		w_ = _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(w_ = _STD move(pack.x)))
		-> decltype(w_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = _STD move(pack.x);
		MXY::y_ = _STD move(pack.y);
		MXYZ::z_ = _STD move(pack.z);
		w_ = _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector4<T>& vector) & noexcept(noexcept(w_ += vector.w_))
		-> decltype(w_ += vector.w_, _STD declval<LRefEucVector>()) {
		MX::x_ += vector.x_;
		MXY::y_ += vector.y_;
		MXYZ::z_ += vector.z_;
		w_ += vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector4<T>&& vector) & noexcept(noexcept(w_ += _STD move(vector.w_)))
		-> decltype(w_ += _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(vector.x_);
		MXY::y_ += _STD move(vector.y_);
		MXYZ::z_ += _STD move(vector.z_);
		w_ += _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(w_ += _STD move(pack.x)))
		-> decltype(w_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ += _STD move(pack.x);
		MXY::y_ += _STD move(pack.y);
		MXYZ::z_ += _STD move(pack.z);
		w_ += _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanRecVector4<T>& vector) & noexcept(noexcept(w_ = w_ + vector.w_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, w_ = w_ + vector.w_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + vector.x_;
		MXY::y_ = MXY::y_ + vector.y_;
		MXYZ::z_ = MXYZ::z_ + vector.z_;
		w_ = w_ + vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanRecVector4<T>&& vector) & noexcept(noexcept(w_ = w_ + _STD move(vector.w_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, w_ = w_ + _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(vector.x_);
		MXY::y_ = MXY::y_ + _STD move(vector.y_);
		MXYZ::z_ = MXYZ::z_ + _STD move(vector.z_);
		w_ = w_ + _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(w_ = w_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, w_ = w_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ + _STD move(pack.x);
		MXY::y_ = MXY::y_ + _STD move(pack.y);
		MXYZ::z_ = MXYZ::z_ + _STD move(pack.z);
		w_ = w_ + _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector4<T>& vector) & noexcept(noexcept(w_ -= vector.w_))
		-> decltype(w_ -= vector.w_, _STD declval<LRefEucVector>()) {
		MX::x_ -= vector.x_;
		MXY::y_ -= vector.y_;
		MXYZ::z_ -= vector.z_;
		w_ -= vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector4<T>&& vector) & noexcept(noexcept(w_ -= _STD move(vector.w_)))
		-> decltype(w_ -= _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(vector.x_);
		MXY::y_ -= _STD move(vector.y_);
		MXYZ::z_ -= _STD move(vector.z_);
		w_ -= _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(w_ -= _STD move(pack.x)))
		-> decltype(w_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ -= _STD move(pack.x);
		MXY::y_ -= _STD move(pack.y);
		MXYZ::z_ -= _STD move(pack.z);
		w_ -= _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanRecVector4<T>& vector) & noexcept(noexcept(w_ = w_ - vector.w_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, w_ = w_ - vector.w_, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - vector.x_;
		MXY::y_ = MXY::y_ - vector.y_;
		MXYZ::z_ = MXYZ::z_ - vector.z_;
		w_ = w_ - vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanRecVector4<T>&& vector) & noexcept(noexcept(w_ = w_ - _STD move(vector.w_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, w_ = w_ - _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(vector.x_);
		MXY::y_ = MXY::y_ - _STD move(vector.y_);
		MXYZ::z_ = MXYZ::z_ - _STD move(vector.z_);
		w_ = w_ - _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(w_ = w_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, w_ = w_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ - _STD move(pack.x);
		MXY::y_ = MXY::y_ - _STD move(pack.y);
		MXYZ::z_ = MXYZ::z_ - _STD move(pack.z);
		w_ = w_ - _STD move(pack.w);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(w_ *= scl))
		-> decltype(w_ *= scl, _STD declval<LRefEucVector>()) {
		MX::x_ *= scl;
		MXY::y_ *= scl;
		MXYZ::z_ *= scl;
		w_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(w_ = w_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, w_ = w_ * scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ * scl;
		MXY::y_ = MXY::y_ * scl;
		MXYZ::z_ = MXYZ::z_ * scl;
		w_ = w_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(w_ /= scl))
		-> decltype(w_ /= scl, _STD declval<LRefEucVector>()) {
		MX::x_ /= scl;
		MXY::y_ /= scl;
		MXYZ::z_ /= scl;
		w_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(w_ = w_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, w_ = w_ / scl, _STD declval<LRefEucVector>()) {
		MX::x_ = MX::x_ / scl;
		MXY::y_ = MXY::y_ / scl;
		MXYZ::z_ = MXYZ::z_ / scl;
		w_ = w_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
#ifdef _MSC_VER
#pragma region 
#endif
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return MX::x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return MXY::y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return  MXY::y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType z() noexcept { return MXYZ::z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType z() const noexcept { return MXYZ::z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType w() noexcept { return w_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType w() const noexcept { return w_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr MXY& xy() noexcept { return *this; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr const MXY& xy() const noexcept { return *this; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr MXYZ& xyz() noexcept { return *this; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr const MXYZ& xyz() const noexcept { return *this; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xx() const noexcept { return { MX::x_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xz() const noexcept { return { MX::x_,MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xw() const noexcept { return { MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yx() const noexcept { return { MXY::y_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yy() const noexcept { return { MXY::y_,MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yz() const noexcept { return { MXY::y_,MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yw() const noexcept { return { MXY::y_,w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zx() const noexcept { return { MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zy() const noexcept { return { MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zz() const noexcept { return { MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zw() const noexcept { return { MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wx() const noexcept { return { w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wy() const noexcept { return { w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wz() const noexcept { return { w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> ww() const noexcept { return { w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxx() const noexcept { return { MX::x_,MX::x_,MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxy() const noexcept { return { MX::x_,MX::x_,MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxz() const noexcept { return { MX::x_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxw() const noexcept { return { MX::x_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyx() const noexcept { return { MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyy() const noexcept { return { MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyw() const noexcept { return { MX::x_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzx() const noexcept { return { MX::x_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzy() const noexcept { return { MX::x_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzz() const noexcept { return { MX::x_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzw() const noexcept { return { MX::x_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwx() const noexcept { return { MX::x_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwy() const noexcept { return { MX::x_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwz() const noexcept { return { MX::x_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xww() const noexcept { return { MX::x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxx() const noexcept { return { MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxy() const noexcept { return { MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxz() const noexcept { return { MXY::y_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxw() const noexcept { return { MXY::y_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyx() const noexcept { return { MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyy() const noexcept { return { MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyz() const noexcept { return { MXY::y_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyw() const noexcept { return { MXY::y_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzx() const noexcept { return { MXY::y_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzy() const noexcept { return { MXY::y_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzz() const noexcept { return { MXY::y_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzw() const noexcept { return { MXY::y_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywx() const noexcept { return { MXY::y_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywy() const noexcept { return { MXY::y_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywz() const noexcept { return { MXY::y_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yww() const noexcept { return { MXY::y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxx() const noexcept { return { MXYZ::z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxy() const noexcept { return { MXYZ::z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxz() const noexcept { return { MXYZ::z_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxw() const noexcept { return { MXYZ::z_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyx() const noexcept { return { MXYZ::z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyy() const noexcept { return { MXYZ::z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyz() const noexcept { return { MXYZ::z_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyw() const noexcept { return { MXYZ::z_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzx() const noexcept { return { MXYZ::z_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzy() const noexcept { return { MXYZ::z_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzz() const noexcept { return { MXYZ::z_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzw() const noexcept { return { MXYZ::z_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwx() const noexcept { return { MXYZ::z_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwy() const noexcept { return { MXYZ::z_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwz() const noexcept { return { MXYZ::z_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zww() const noexcept { return { MXYZ::z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxx() const noexcept { return { w_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxy() const noexcept { return { w_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxz() const noexcept { return { w_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxw() const noexcept { return { w_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyx() const noexcept { return { w_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyy() const noexcept { return { w_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyz() const noexcept { return { w_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyw() const noexcept { return { w_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzx() const noexcept { return { w_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzy() const noexcept { return { w_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzz() const noexcept { return { w_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzw() const noexcept { return { w_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwx() const noexcept { return { w_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwy() const noexcept { return { w_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwz() const noexcept { return { w_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> www() const noexcept { return { w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxx() const noexcept { return { MX::x_, MX::x_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxy() const noexcept { return { MX::x_, MX::x_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxz() const noexcept { return { MX::x_, MX::x_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxw() const noexcept { return { MX::x_, MX::x_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyx() const noexcept { return { MX::x_, MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyy() const noexcept { return { MX::x_, MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyz() const noexcept { return { MX::x_, MX::x_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyw() const noexcept { return { MX::x_, MX::x_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzx() const noexcept { return { MX::x_, MX::x_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzy() const noexcept { return { MX::x_, MX::x_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzz() const noexcept { return { MX::x_, MX::x_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzw() const noexcept { return { MX::x_, MX::x_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwx() const noexcept { return { MX::x_, MX::x_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwy() const noexcept { return { MX::x_, MX::x_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwz() const noexcept { return { MX::x_, MX::x_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxww() const noexcept { return { MX::x_, MX::x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxx() const noexcept { return { MX::x_, MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxy() const noexcept { return { MX::x_, MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxz() const noexcept { return { MX::x_, MXY::y_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxw() const noexcept { return { MX::x_, MXY::y_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyx() const noexcept { return { MX::x_, MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyy() const noexcept { return { MX::x_, MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyz() const noexcept { return { MX::x_, MXY::y_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyw() const noexcept { return { MX::x_, MXY::y_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzx() const noexcept { return { MX::x_, MXY::y_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzy() const noexcept { return { MX::x_, MXY::y_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzz() const noexcept { return { MX::x_, MXY::y_, MXYZ::z_, MXYZ::z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywx() const noexcept { return { MX::x_, MXY::y_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywy() const noexcept { return { MX::x_, MXY::y_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywz() const noexcept { return { MX::x_, MXY::y_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyww() const noexcept { return { MX::x_, MXY::y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxx() const noexcept { return { MX::x_, MXYZ::z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxy() const noexcept { return { MX::x_, MXYZ::z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxz() const noexcept { return { MX::x_, MXYZ::z_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxw() const noexcept { return { MX::x_, MXYZ::z_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyx() const noexcept { return { MX::x_, MXYZ::z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyy() const noexcept { return { MX::x_, MXYZ::z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyz() const noexcept { return { MX::x_, MXYZ::z_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyw() const noexcept { return { MX::x_, MXYZ::z_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzx() const noexcept { return { MX::x_, MXYZ::z_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzy() const noexcept { return { MX::x_, MXYZ::z_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzz() const noexcept { return { MX::x_, MXYZ::z_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzw() const noexcept { return { MX::x_, MXYZ::z_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwx() const noexcept { return { MX::x_, MXYZ::z_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwy() const noexcept { return { MX::x_, MXYZ::z_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwz() const noexcept { return { MX::x_, MXYZ::z_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzww() const noexcept { return { MX::x_, MXYZ::z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxx() const noexcept { return { MX::x_, w_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxy() const noexcept { return { MX::x_, w_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxz() const noexcept { return { MX::x_, w_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxw() const noexcept { return { MX::x_, w_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyx() const noexcept { return { MX::x_, w_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyy() const noexcept { return { MX::x_, w_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyz() const noexcept { return { MX::x_, w_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyw() const noexcept { return { MX::x_, w_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzx() const noexcept { return { MX::x_, w_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzy() const noexcept { return { MX::x_, w_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzz() const noexcept { return { MX::x_, w_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzw() const noexcept { return { MX::x_, w_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwx() const noexcept { return { MX::x_, w_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwy() const noexcept { return { MX::x_, w_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwz() const noexcept { return { MX::x_, w_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwww() const noexcept { return { MX::x_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxx() const noexcept { return { MXY::y_, MX::x_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxy() const noexcept { return { MXY::y_, MX::x_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxz() const noexcept { return { MXY::y_, MX::x_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxw() const noexcept { return { MXY::y_, MX::x_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyx() const noexcept { return { MXY::y_, MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyy() const noexcept { return { MXY::y_, MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyz() const noexcept { return { MXY::y_, MX::x_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyw() const noexcept { return { MXY::y_, MX::x_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzx() const noexcept { return { MXY::y_, MX::x_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzy() const noexcept { return { MXY::y_, MX::x_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzz() const noexcept { return { MXY::y_, MX::x_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzw() const noexcept { return { MXY::y_, MX::x_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwx() const noexcept { return { MXY::y_, MX::x_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwy() const noexcept { return { MXY::y_, MX::x_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwz() const noexcept { return { MXY::y_, MX::x_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxww() const noexcept { return { MXY::y_, MX::x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxx() const noexcept { return { MXY::y_, MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxy() const noexcept { return { MXY::y_, MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxz() const noexcept { return { MXY::y_, MXY::y_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxw() const noexcept { return { MXY::y_, MXY::y_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyx() const noexcept { return { MXY::y_, MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyy() const noexcept { return { MXY::y_, MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyz() const noexcept { return { MXY::y_, MXY::y_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyw() const noexcept { return { MXY::y_, MXY::y_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzx() const noexcept { return { MXY::y_, MXY::y_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzy() const noexcept { return { MXY::y_, MXY::y_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzz() const noexcept { return { MXY::y_, MXY::y_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzw() const noexcept { return { MXY::y_, MXY::y_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywx() const noexcept { return { MXY::y_, MXY::y_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywy() const noexcept { return { MXY::y_, MXY::y_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywz() const noexcept { return { MXY::y_, MXY::y_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyww() const noexcept { return { MXY::y_, MXY::y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxx() const noexcept { return { MXY::y_, MXYZ::z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxy() const noexcept { return { MXY::y_, MXYZ::z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxz() const noexcept { return { MXY::y_, MXYZ::z_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxw() const noexcept { return { MXY::y_, MXYZ::z_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyx() const noexcept { return { MXY::y_, MXYZ::z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyy() const noexcept { return { MXY::y_, MXYZ::z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyz() const noexcept { return { MXY::y_, MXYZ::z_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyw() const noexcept { return { MXY::y_, MXYZ::z_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzx() const noexcept { return { MXY::y_, MXYZ::z_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzy() const noexcept { return { MXY::y_, MXYZ::z_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzz() const noexcept { return { MXY::y_, MXYZ::z_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzw() const noexcept { return { MXY::y_, MXYZ::z_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwx() const noexcept { return { MXY::y_, MXYZ::z_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwy() const noexcept { return { MXY::y_, MXYZ::z_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwz() const noexcept { return { MXY::y_, MXYZ::z_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzww() const noexcept { return { MXY::y_, MXYZ::z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxx() const noexcept { return { MXY::y_, w_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxy() const noexcept { return { MXY::y_, w_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxz() const noexcept { return { MXY::y_, w_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxw() const noexcept { return { MXY::y_, w_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyx() const noexcept { return { MXY::y_, w_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyy() const noexcept { return { MXY::y_, w_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyz() const noexcept { return { MXY::y_, w_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyw() const noexcept { return { MXY::y_, w_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzx() const noexcept { return { MXY::y_, w_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzy() const noexcept { return { MXY::y_, w_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzz() const noexcept { return { MXY::y_, w_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzw() const noexcept { return { MXY::y_, w_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwx() const noexcept { return { MXY::y_, w_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwy() const noexcept { return { MXY::y_, w_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwz() const noexcept { return { MXY::y_, w_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywww() const noexcept { return { MXY::y_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxx() const noexcept { return { MXYZ::z_, MX::x_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxy() const noexcept { return { MXYZ::z_, MX::x_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxz() const noexcept { return { MXYZ::z_, MX::x_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxw() const noexcept { return { MXYZ::z_, MX::x_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyx() const noexcept { return { MXYZ::z_, MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyy() const noexcept { return { MXYZ::z_, MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyz() const noexcept { return { MXYZ::z_, MX::x_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyw() const noexcept { return { MXYZ::z_, MX::x_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzx() const noexcept { return { MXYZ::z_, MX::x_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzy() const noexcept { return { MXYZ::z_, MX::x_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzz() const noexcept { return { MXYZ::z_, MX::x_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzw() const noexcept { return { MXYZ::z_, MX::x_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwx() const noexcept { return { MXYZ::z_, MX::x_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwy() const noexcept { return { MXYZ::z_, MX::x_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwz() const noexcept { return { MXYZ::z_, MX::x_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxww() const noexcept { return { MXYZ::z_, MX::x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxx() const noexcept { return { MXYZ::z_, MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxy() const noexcept { return { MXYZ::z_, MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxz() const noexcept { return { MXYZ::z_, MXY::y_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxw() const noexcept { return { MXYZ::z_, MXY::y_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyx() const noexcept { return { MXYZ::z_, MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyy() const noexcept { return { MXYZ::z_, MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyz() const noexcept { return { MXYZ::z_, MXY::y_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyw() const noexcept { return { MXYZ::z_, MXY::y_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzx() const noexcept { return { MXYZ::z_, MXY::y_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzy() const noexcept { return { MXYZ::z_, MXY::y_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzz() const noexcept { return { MXYZ::z_, MXY::y_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzw() const noexcept { return { MXYZ::z_, MXY::y_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywx() const noexcept { return { MXYZ::z_, MXY::y_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywy() const noexcept { return { MXYZ::z_, MXY::y_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywz() const noexcept { return { MXYZ::z_, MXY::y_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyww() const noexcept { return { MXYZ::z_, MXY::y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxx() const noexcept { return { MXYZ::z_, MXYZ::z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxy() const noexcept { return { MXYZ::z_, MXYZ::z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxz() const noexcept { return { MXYZ::z_, MXYZ::z_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxw() const noexcept { return { MXYZ::z_, MXYZ::z_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyx() const noexcept { return { MXYZ::z_, MXYZ::z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyy() const noexcept { return { MXYZ::z_, MXYZ::z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyz() const noexcept { return { MXYZ::z_, MXYZ::z_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyw() const noexcept { return { MXYZ::z_, MXYZ::z_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzx() const noexcept { return { MXYZ::z_, MXYZ::z_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzy() const noexcept { return { MXYZ::z_, MXYZ::z_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzz() const noexcept { return { MXYZ::z_, MXYZ::z_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzw() const noexcept { return { MXYZ::z_, MXYZ::z_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwx() const noexcept { return { MXYZ::z_, MXYZ::z_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwy() const noexcept { return { MXYZ::z_, MXYZ::z_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwz() const noexcept { return { MXYZ::z_, MXYZ::z_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzww() const noexcept { return { MXYZ::z_, MXYZ::z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxx() const noexcept { return { MXYZ::z_, w_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxy() const noexcept { return { MXYZ::z_, w_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxz() const noexcept { return { MXYZ::z_, w_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxw() const noexcept { return { MXYZ::z_, w_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyx() const noexcept { return { MXYZ::z_, w_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyy() const noexcept { return { MXYZ::z_, w_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyz() const noexcept { return { MXYZ::z_, w_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyw() const noexcept { return { MXYZ::z_, w_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzx() const noexcept { return { MXYZ::z_, w_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzy() const noexcept { return { MXYZ::z_, w_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzz() const noexcept { return { MXYZ::z_, w_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzw() const noexcept { return { MXYZ::z_, w_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwx() const noexcept { return { MXYZ::z_, w_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwy() const noexcept { return { MXYZ::z_, w_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwz() const noexcept { return { MXYZ::z_, w_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwww() const noexcept { return { MXYZ::z_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxx() const noexcept { return { w_, MX::x_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxy() const noexcept { return { w_, MX::x_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxz() const noexcept { return { w_, MX::x_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxw() const noexcept { return { w_, MX::x_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyx() const noexcept { return { w_, MX::x_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyy() const noexcept { return { w_, MX::x_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyz() const noexcept { return { w_, MX::x_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyw() const noexcept { return { w_, MX::x_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzx() const noexcept { return { w_, MX::x_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzy() const noexcept { return { w_, MX::x_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzz() const noexcept { return { w_, MX::x_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzw() const noexcept { return { w_, MX::x_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwx() const noexcept { return { w_, MX::x_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwy() const noexcept { return { w_, MX::x_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwz() const noexcept { return { w_, MX::x_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxww() const noexcept { return { w_, MX::x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxx() const noexcept { return { w_, MXY::y_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxy() const noexcept { return { w_, MXY::y_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxz() const noexcept { return { w_, MXY::y_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxw() const noexcept { return { w_, MXY::y_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyx() const noexcept { return { w_, MXY::y_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyy() const noexcept { return { w_, MXY::y_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyz() const noexcept { return { w_, MXY::y_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyw() const noexcept { return { w_, MXY::y_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzx() const noexcept { return { w_, MXY::y_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzy() const noexcept { return { w_, MXY::y_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzz() const noexcept { return { w_, MXY::y_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzw() const noexcept { return { w_, MXY::y_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywx() const noexcept { return { w_, MXY::y_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywy() const noexcept { return { w_, MXY::y_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywz() const noexcept { return { w_, MXY::y_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyww() const noexcept { return { w_, MXY::y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxx() const noexcept { return { w_, MXYZ::z_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxy() const noexcept { return { w_, MXYZ::z_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxz() const noexcept { return { w_, MXYZ::z_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxw() const noexcept { return { w_, MXYZ::z_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyx() const noexcept { return { w_, MXYZ::z_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyy() const noexcept { return { w_, MXYZ::z_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyz() const noexcept { return { w_, MXYZ::z_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyw() const noexcept { return { w_, MXYZ::z_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzx() const noexcept { return { w_, MXYZ::z_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzy() const noexcept { return { w_, MXYZ::z_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzz() const noexcept { return { w_, MXYZ::z_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzw() const noexcept { return { w_, MXYZ::z_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwx() const noexcept { return { w_, MXYZ::z_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwy() const noexcept { return { w_, MXYZ::z_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwz() const noexcept { return { w_, MXYZ::z_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzww() const noexcept { return { w_, MXYZ::z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxx() const noexcept { return { w_, w_, MX::x_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxy() const noexcept { return { w_, w_, MX::x_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxz() const noexcept { return { w_, w_, MX::x_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxw() const noexcept { return { w_, w_, MX::x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyx() const noexcept { return { w_, w_, MXY::y_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyy() const noexcept { return { w_, w_, MXY::y_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyz() const noexcept { return { w_, w_, MXY::y_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyw() const noexcept { return { w_, w_, MXY::y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzx() const noexcept { return { w_, w_, MXYZ::z_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzy() const noexcept { return { w_, w_, MXYZ::z_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzz() const noexcept { return { w_, w_, MXYZ::z_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzw() const noexcept { return { w_, w_, MXYZ::z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwx() const noexcept { return { w_, w_, w_, MX::x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwy() const noexcept { return { w_, w_, w_, MXY::y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwz() const noexcept { return { w_, w_, w_, MXYZ::z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwww() const noexcept { return { w_, w_, w_, w_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { MX::x_,MXY::y_, MXYZ::z_, w_ }; }
#ifdef _MSC_VER
#pragma endregion
#endif
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		MX::x_ = 0;
		MXY::y_ = 0;
		MXYZ::z_ = 0;
		w_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y, class Z, class W>
	EUCVECTORINLINE auto set(X&& x, Y&& y, Z&& z, W&& w) noexcept(noexcept(MXYZ::set(_STD forward<X>(x), _STD forward<Y>(y), _STD forward<Z>(z))) && noexcept(w_ = _STD forward<W>(w)))
		-> meta::is_type_t<decltype(MXYZ::set(_STD forward<X>(x), _STD forward<Y>(y), _STD forward<Z>(z)), w_ = _STD forward<W>(w)), void> {
		MX::x_ = _STD forward<X>(x);
		MXY::y_ = _STD forward<Y>(y);
		MXYZ::z_ = _STD forward<Z>(z);
		w_ = _STD forward<W>(w);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanRecVector4<T>& vector) const noexcept(noexcept(MXYZ::dot(vector) + decltype(w_ * vector.w_)(w_ * vector.w_)))
		-> decltype(MXYZ::dot(vector) + decltype(w_ * vector.w_)(w_ * vector.w_)) {
		return (MX::x_ * vector.x_) + (MXY::y_ * vector.y_) + (MXYZ::z_ * vector.z_) + (w_ * vector.w_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanRecVector4<T>&& vector) const noexcept(noexcept(MXYZ::dot(_STD move(vector)) + decltype(w_ * _STD move(vector.w_))(w_ * _STD move(vector.w_))))
		-> decltype(MXYZ::dot(_STD move(vector)) + decltype(w_ * _STD move(vector.w_))(w_ * _STD move(vector.w_))) {
		return (MX::x_ * _STD move(vector.x_)) + (MXY::y_ * _STD move(vector.y_)) + (MXYZ::z_ * _STD move(vector.z_)) + (w_ * _STD move(vector.w_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(MXYZ::dot(_STD move(pack)) + decltype(w_ * _STD move(pack.x))(w_ * _STD move(pack.x))))
		-> decltype(MXYZ::dot(_STD move(pack)) + decltype(w_ * _STD move(pack.x))(w_ * _STD move(pack.x))) {
		return (MX::x_ * _STD move(pack.x)) + (MXY::y_ * _STD move(pack.y)) + (MXYZ::z_ * _STD move(pack.z)) + (w_ * _STD move(pack.w));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return (MX::x_ * MX::x_) + (MXY::y_ * MXY::y_) + (MXYZ::z_ * MXYZ::z_) + (w_ * w_);
	}
	/*
		@brief

			Calculate the norm.

			out = ROOT((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt((MX::x_ * MX::x_) + (MXY::y_ * MXY::y_) + (MXYZ::z_ * MXYZ::z_) + (w_ * w_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{w_ / eucnorm<T>()}))
		->decltype(Packer<T>{w_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { MX::x_ / norm, MXY::y_ / norm, MXYZ::z_ / norm, w_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}

};

/*
	Completion Vector.
*/

/*
	D2.
*/
template<class E = float>
struct EuclideanCmplVector2 final
	: private meta::evd_euc_vec {
protected:

	static constexpr size_t EucD = 2;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using ElemType = E;
	using EucVector = EuclideanCmplVector2<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_2<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanCmplVector2;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

public:

	ElemType x_, y_;


	/*
		Constructors.
	*/
	EuclideanCmplVector2() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: x_()
		, y_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanCmplVector2(const EuclideanCmplVector2& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: x_(vector.x_)
		, y_(vector.y_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanCmplVector2(EuclideanCmplVector2&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: x_(_STD move(vector.x_))
		, y_(_STD move(vector.y_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanCmplVector2(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: x_(_STD move(pack.x))
		, y_(_STD move(pack.y))
	{}

	template<class X, class Y, meta::if_t<meta::is_constructible_anynum_param_v<ElemType, X, Y>> = 0>
	EuclideanCmplVector2(X&& x, Y&& y) noexcept(noexcept(MX(_STD forward<X>(x))) && _STD is_nothrow_constructible_v<ElemType, Y>)
		: x_(_STD forward<X>(x))
		, y_(_STD forward<Y>(y))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanCmplVector2<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + vector.y_)>>{y_ + vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(y_ + vector.y_)>>{y_ + vector.y_}) {
		return { x_ + vector.x_, y_ + vector.y_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanCmplVector2<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + _STD move(vector.y_))>>{y_ + _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(y_ + _STD move(vector.y_))>>{y_ + _STD move(vector.y_)}) {
		return { x_ + _STD move(vector.x_), y_ + _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ + _STD move(pack.x))>>{y_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(y_ + _STD move(pack.x))>>{y_ + _STD move(pack.x)}) {
		return { x_ + _STD move(pack.x), y_ + _STD move(pack.y) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.y_)>>{_STD move(pack.x) + vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.y_)>>{_STD move(pack.x) + vector.y_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.y_))>>{_STD move(pack.x) + _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.y_))>>{_STD move(pack.x) + _STD move(vector.y_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_),_STD move(pack.y) + _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanCmplVector2<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - vector.y_)>>{y_ - vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(y_ - vector.y_)>>{y_ - vector.y_}) {
		return { x_ - vector.x_, y_ - vector.y_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanCmplVector2<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - _STD move(vector.y_))>>{y_ - _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(y_ - _STD move(vector.y_))>>{y_ - _STD move(vector.y_)}) {
		return { x_ - _STD move(vector.x_), y_ - _STD move(vector.y_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ - _STD move(pack.x))>>{y_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(y_ - _STD move(pack.x))>>{y_ - _STD move(pack.x)}) {
		return { x_ - _STD move(pack.x), y_ - _STD move(pack.y) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.y_)>>{_STD move(pack.x) - vector.y_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.y_)>>{_STD move(pack.x) - vector.y_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.y_))>>{_STD move(pack.x) - _STD move(vector.y_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.y_))>>{_STD move(pack.x) - _STD move(vector.y_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_),_STD move(pack.y) - _STD move(vector.y_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_* scl)>>{y_* scl}))
		->decltype(Packer<meta::no_ref<decltype(y_* scl)>>{y_* scl}) {
		return { x_ * scl, y_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.y_* scl)>>{right.y_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.y_* scl)>>{right.y_* scl}) {
		return { right.x_ * scl, right.y_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.y_)* scl)>>{_STD move(right.y_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.y_)* scl)>>{_STD move(right.y_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(y_ / scl)>>{y_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(y_ / scl)>>{y_ / scl}) {
		return { x_ / scl, y_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanCmplVector2<T>& right) const noexcept(noexcept(bool(y_ == right.y_)))
		-> meta::no_ref<decltype(bool(y_ == right.y_))> {
		return (x_ == right.x_) && (y_ == right.y_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanCmplVector2<T>&& right) const noexcept(noexcept(bool(y_ == _STD move(right.y_))))
		-> meta::no_ref<decltype(bool(y_ == _STD move(right.y_)))> {
		return (x_ == _STD move(right.x_)) && (y_ == _STD move(right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(y_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(y_ == _STD move(right.x)))> {
		return (x_ == _STD move(right.x)) && (y_ == _STD move(right.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.y_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.y_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.y_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.y_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector2<T>& right) const noexcept(noexcept(bool(y_ != right.y_)))
		-> meta::no_ref<decltype(bool(y_ != right.y_))> {
		return (x_ != right.x_) || (y_ != right.y_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector2<T>&& right) const noexcept(noexcept(bool(y_ != _STD move(right.y_))))
		-> meta::no_ref<decltype(bool(y_ != _STD move(right.y_)))> {
		return (x_ != _STD move(right.x_)) || (y_ != _STD move(right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(y_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(y_ != _STD move(right.x)), _STD declval<bool>())> {
		return (x_ != _STD move(right.x)) || (y_ != _STD move(right.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.y_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.y_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.y_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.y_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector2<T>& right) const noexcept(noexcept(bool(!(y_ == right.y_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(y_ == right.y_)), _STD declval<bool>())> {
		return (!(x_ == right.x_)) || (!(y_ == right.y_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector2<T>&& right) const noexcept(noexcept(bool(!(y_ == _STD move(right.y_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(y_ == _STD move(right.y_))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x_))) || (!(y_ == _STD move(right.y_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(y_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(y_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x))) || (!(y_ == _STD move(right.y)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.y_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.y_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.y_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanCmplVector2<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanCmplVector2<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanCmplVector2<T>() * -1))
		->decltype(EuclideanCmplVector2<T>() * -1) {
		return { x_ * -1, y_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanCmplVector2<T>& vector) & noexcept(noexcept(y_ = vector.y_))
		-> decltype(y_ = vector.y_, _STD declval<LRefEucVector>()) {
		x_ = vector.x_;
		y_ = vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanCmplVector2<T>&& vector) & noexcept(noexcept(y_ = _STD move(vector.y_)))
		-> decltype(y_ = _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		x_ = _STD move(vector.x_);
		y_ = _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(y_ = _STD move(pack.x)))
		-> decltype(y_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = _STD move(pack.x);
		y_ = _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector2<T>& vector) & noexcept(noexcept(y_ += vector.y_))
		-> decltype(y_ += vector.y_, _STD declval<LRefEucVector>()) {
		x_ += vector.x_;
		y_ += vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector2<T>&& vector) & noexcept(noexcept(y_ += _STD move(vector.y_)))
		-> decltype(y_ += _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		x_ += _STD move(vector.x_);
		y_ += _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(y_ += _STD move(pack.x)))
		-> decltype(y_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ += _STD move(pack.x);
		y_ += _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector2<T>& vector) & noexcept(noexcept(y_ = y_ + vector.y_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, y_ = y_ + vector.y_, _STD declval<LRefEucVector>()) {
		x_ = x_ + vector.x_;
		y_ = y_ + vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector2<T>&& vector) & noexcept(noexcept(y_ = y_ + _STD move(vector.y_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, y_ = y_ + _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(vector.x_);
		y_ = y_ + _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(y_ = y_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, y_ = y_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(pack.x);
		y_ = y_ + _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector2<T>& vector) & noexcept(noexcept(y_ -= vector.y_))
		-> decltype(y_ -= vector.y_, _STD declval<LRefEucVector>()) {
		x_ -= vector.x_;
		y_ -= vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector2<T>&& vector) & noexcept(noexcept(y_ -= _STD move(vector.y_)))
		-> decltype(y_ -= _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(vector.x_);
		y_ -= _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(y_ -= _STD move(pack.x)))
		-> decltype(y_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(pack.x);
		y_ -= _STD move(pack.y);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector2<T>& vector) & noexcept(noexcept(y_ = y_ - vector.y_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, y_ = y_ - vector.y_, _STD declval<LRefEucVector>()) {
		x_ = x_ - vector.x_;
		y_ = y_ - vector.y_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector2<T>&& vector) & noexcept(noexcept(y_ = y_ - _STD move(vector.y_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, y_ = y_ - _STD move(vector.y_), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(vector.x_);
		y_ = y_ - _STD move(vector.y_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(y_ = y_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, y_ = y_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(pack.x);
		y_ = y_ - _STD move(pack.y);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(y_ *= scl))
		-> decltype(y_ *= scl, _STD declval<LRefEucVector>()) {
		x_ *= scl;
		y_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(y_ = y_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, y_ = y_ * scl, _STD declval<LRefEucVector>()) {
		x_ = x_ * scl;
		y_ = y_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(y_ /= scl))
		-> decltype(y_ /= scl, _STD declval<LRefEucVector>()) {
		x_ /= scl;
		y_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(y_ = y_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, y_ = y_ / scl, _STD declval<LRefEucVector>()) {
		x_ = x_ / scl;
		y_ = y_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xx() const noexcept { return { x_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yx() const noexcept { return { y_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yy() const noexcept { return { y_,y_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { x_,y_ }; }
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		x_ = 0;
		y_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y>
	EUCVECTORINLINE auto set(X&& x, Y&& y) noexcept(noexcept(x_ = _STD forward<X>(x)) && noexcept(y_ = _STD forward<Y>(y)))
		-> meta::is_type_t<decltype(x_ = _STD forward<X>(x), y_ = _STD forward<Y>(y)), void> {
		x_ = _STD forward<X>(x);
		y_ = _STD forward<Y>(y);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanCmplVector2<T>& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(x_* vector.x_ + y_* vector.y_)))
		-> decltype(ElemType(x_ * vector.x_ + y_ * vector.y_), x_ + x_) {
		return ElemType(x_ * vector.x_) + ElemType(y_ * vector.y_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanCmplVector2<T>&& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(x_* _STD move(vector.x_) + y_ * _STD move(vector.y_))))
		-> decltype(ElemType(x_* _STD move(vector.x_) + y_ * _STD move(vector.y_)), x_ + x_) {
		return ElemType(x_ * _STD move(vector.x_)) + ElemType(y_ * _STD move(vector.y_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(x_* _STD move(pack.x) + y_ * _STD move(pack.y))))
		-> decltype(ElemType(x_* _STD move(pack.x) + y_ * _STD move(pack.y)), x_ + x_) {
		return ElemType(x_ * _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return T(x_ * x_) + T(y_ * y_);
	}
	/*
		@brief

			Calculate the norm.

			out = ROOT((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt(T(x_ * x_) + T(y_ * y_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{y_ / eucnorm<T>()}))
		->decltype(Packer<T>{y_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { x_ / norm, y_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}

};

/*
	D3.
*/
template<class E = float>
struct EuclideanCmplVector3 final
	: private meta::evd_euc_vec {
protected:

	static constexpr size_t EucD = 3;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using ElemType = E;
	using EucVector = EuclideanCmplVector3<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_3<R>;
	template<class R>
	using SubPacker2 = detail::ResultPacker_2<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanCmplVector3;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

public:

	ElemType x_, y_, z_;

	/*
		Constructors.
	*/
	EuclideanCmplVector3() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: x_()
		, y_()
		, z_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanCmplVector3(const EuclideanCmplVector3& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: x_(vector.x_)
		, y_(vector.y_)
		, z_(vector.z_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanCmplVector3(EuclideanCmplVector3&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: x_(_STD move(vector.x_))
		, y_(_STD move(vector.y_))
		, z_(_STD move(vector.z_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanCmplVector3(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: x_(_STD move(pack.x_))
		, y_(_STD move(pack.y_))
		, z_(_STD move(pack.z_))
	{}

	template<class X, class Y, class Z, meta::if_t<meta::is_constructible_anynum_param_v<ElemType, X, Y, Z>> = 0>
	EuclideanCmplVector3(X&& x, Y&& y, Z&& z) noexcept(_STD is_nothrow_constructible_v<ElemType, X> && _STD is_nothrow_constructible_v<ElemType, Y> && _STD is_nothrow_constructible_v<ElemType, Z>)
		: x_(_STD forward<X>(x))
		, y_(_STD forward<Y>(y))
		, z_(_STD forward<Z>(z))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanCmplVector3<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + vector.z_)>>{z_ + vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(z_ + vector.z_)>>{z_ + vector.z_}) {
		return { x_ + vector.x_, y_ + vector.y_, z_ + vector.z_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanCmplVector3<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + _STD move(vector.z_))>>{z_ + _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(z_ + _STD move(vector.z_))>>{z_ + _STD move(vector.z_)}) {
		return { x_ + _STD move(vector.x_), y_ + _STD move(vector.y_), z_ + _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ + _STD move(pack.x))>>{z_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(z_ + _STD move(pack.x))>>{z_ + _STD move(pack.x)}) {
		return { x_ + _STD move(pack.x), y_ + _STD move(pack.y), z_ + _STD move(pack.z) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.z_)>>{_STD move(pack.x) + vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.z_)>>{_STD move(pack.x) + vector.z_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ , _STD move(pack.z) + vector.z_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.z_))>>{_STD move(pack.x) + _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.z_))>>{_STD move(pack.x) + _STD move(vector.z_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_), _STD move(pack.y) + _STD move(vector.y_), _STD move(pack.z) + _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanCmplVector3<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - vector.z_)>>{z_ - vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(z_ - vector.z_)>>{z_ - vector.z_}) {
		return { x_ - vector.x_, y_ - vector.y_, z_ - vector.z_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanCmplVector3<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - _STD move(vector.z_))>>{z_ - _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(z_ - _STD move(vector.z_))>>{z_ - _STD move(vector.z_)}) {
		return { x_ - _STD move(vector.x_), y_ - _STD move(vector.y_), z_ - _STD move(vector.z_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ - _STD move(pack.x))>>{z_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(z_ - _STD move(pack.x))>>{z_ - _STD move(pack.x)}) {
		return { x_ - _STD move(pack.x), y_ - _STD move(pack.y), z_ - _STD move(pack.z) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.z_)>>{_STD move(pack.x) - vector.z_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.z_)>>{_STD move(pack.x) - vector.z_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ , _STD move(pack.z) - vector.z_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.z_))>>{_STD move(pack.x) - _STD move(vector.z_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.z_))>>{_STD move(pack.x) - _STD move(vector.z_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_), _STD move(pack.y) - _STD move(vector.y_), _STD move(pack.z) - _STD move(vector.z_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_* scl)>>{z_* scl}))
		->decltype(Packer<meta::no_ref<decltype(z_* scl)>>{z_* scl}) {
		return { x_ * scl, y_ * scl, z_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.z_* scl)>>{right.z_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.z_* scl)>>{right.z_* scl}) {
		return { right.x_ * scl,  right.y_ * scl, right.z_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.z_)* scl)>>{_STD move(right.z_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.z_)* scl)>>{_STD move(right.z_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl, _STD move(right.z_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(z_ / scl)>>{z_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(z_ / scl)>>{z_ / scl}) {
		return { x_ / scl, y_ / scl, z_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanCmplVector3<T>& right) const noexcept(noexcept(bool(z_ == right.z_)))
		-> meta::no_ref<decltype(bool(z_ == right.z_))> {
		return (x_ == right.x_) && (y_ == right.y_) && (z_ == right.z_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanCmplVector3<T>&& right) const noexcept(noexcept(bool(z_ == _STD move(right.z_))))
		-> meta::no_ref<decltype(bool(z_ == _STD move(right.z_)))> {
		return (x_ == _STD move(right.x_)) && (y_ == _STD move(right.y_)) && (z_ == _STD move(right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(z_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(z_ == _STD move(right.x)))> {
		return (x_ == _STD move(right.x)) && (y_ == _STD move(right.y)) && (z_ == _STD move(right.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.z_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.z_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y)) && (right.z_ == _STD move(left.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.z_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.z_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y)) && (_STD move(right.z_) == _STD move(left.z));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector3<T>& right) const noexcept(noexcept(bool(z_ != right.z_)))
		-> meta::no_ref<decltype(bool(z_ != right.z_))> {
		return (x_ != right.x_) || (y_ != right.y_) || (z_ != right.z_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector3<T>&& right) const noexcept(noexcept(bool(z_ != _STD move(right.z_))))
		-> meta::no_ref<decltype(bool(z_ != _STD move(right.z_)))> {
		return (x_ != _STD move(right.x_)) || (y_ != _STD move(right.y_)) || (z_ != _STD move(right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(z_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(z_ != _STD move(right.x)), _STD declval<bool>())> {
		return (x_ != _STD move(right.x)) || (y_ != _STD move(right.y)) || (z_ != _STD move(right.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.z_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.z_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y)) || (right.z_ != _STD move(left.z));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.z_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.z_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y)) || (_STD move(right.z_) != _STD move(left.z));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector3<T>& right) const noexcept(noexcept(bool(!(z_ == right.z_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(z_ == right.z_)), _STD declval<bool>())> {
		return (!(x_ == right.x_)) || (!(y_ == right.y_)) || (!(z_ == right.z_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector3<T>&& right) const noexcept(noexcept(bool(!(z_ == _STD move(right.z_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(z_ == _STD move(right.z_))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x_))) || (!(y_ == _STD move(right.y_))) || (!(z_ == _STD move(right.z_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(z_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(z_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x))) || (!(y_ == _STD move(right.y))) || (!(z_ == _STD move(right.z)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.z_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.z_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y))) || (!(right.z_ == _STD move(left.z)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.z_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y))) || (!(_STD move(right.z_) == _STD move(left.z)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanCmplVector3<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanCmplVector3<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanCmplVector3<T>() * -1))
		->decltype(EuclideanCmplVector3<T>() * -1) {
		return { x_ * -1, y_ * -1, z_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanCmplVector3<T>& vector) & noexcept(noexcept(z_ = vector.z_))
		-> decltype(z_ = vector.z_, _STD declval<LRefEucVector>()) {
		x_ = vector.x_;
		y_ = vector.y_;
		z_ = vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanCmplVector3<T>&& vector) & noexcept(noexcept(z_ = _STD move(vector.z_)))
		-> decltype(z_ = _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		x_ = _STD move(vector.x_);
		y_ = _STD move(vector.y_);
		z_ = _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(z_ = _STD move(pack.x)))
		-> decltype(z_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = _STD move(pack.x);
		y_ = _STD move(pack.y);
		z_ = _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector3<T>& vector) & noexcept(noexcept(z_ += vector.z_))
		-> decltype(z_ += vector.z_, _STD declval<LRefEucVector>()) {
		x_ += vector.x_;
		y_ += vector.y_;
		z_ += vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector3<T>&& vector) & noexcept(noexcept(z_ += _STD move(vector.z_)))
		-> decltype(z_ += _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		x_ += _STD move(vector.x_);
		y_ += _STD move(vector.y_);
		z_ += _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(z_ += _STD move(pack.x)))
		-> decltype(z_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ += _STD move(pack.x);
		y_ += _STD move(pack.y);
		z_ += _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector3<T>& vector) & noexcept(noexcept(z_ = z_ + vector.z_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, z_ = z_ + vector.z_, _STD declval<LRefEucVector>()) {
		x_ = x_ + vector.x_;
		y_ = y_ + vector.y_;
		z_ = z_ + vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector3<T>&& vector) & noexcept(noexcept(z_ = z_ + _STD move(vector.z_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, z_ = z_ + _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(vector.x_);
		y_ = y_ + _STD move(vector.y_);
		z_ = z_ + _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(z_ = z_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, z_ = z_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(pack.x);
		y_ = y_ + _STD move(pack.y);
		z_ = z_ + _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector3<T>& vector) & noexcept(noexcept(z_ -= vector.z_))
		-> decltype(z_ -= vector.z_, _STD declval<LRefEucVector>()) {
		x_ -= vector.x_;
		y_ -= vector.y_;
		z_ -= vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector3<T>&& vector) & noexcept(noexcept(z_ -= _STD move(vector.z_)))
		-> decltype(z_ -= _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(vector.x_);
		y_ -= _STD move(vector.y_);
		z_ -= _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(z_ -= _STD move(pack.x)))
		-> decltype(z_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(pack.x);
		y_ -= _STD move(pack.y);
		z_ -= _STD move(pack.z);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector3<T>& vector) & noexcept(noexcept(z_ = z_ - vector.z_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, z_ = z_ - vector.z_, _STD declval<LRefEucVector>()) {
		x_ = x_ - vector.x_;
		y_ = y_ - vector.y_;
		z_ = z_ - vector.z_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector3<T>&& vector) & noexcept(noexcept(z_ = z_ - _STD move(vector.z_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, z_ = z_ - _STD move(vector.z_), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(vector.x_);
		y_ = y_ - _STD move(vector.y_);
		z_ = z_ - _STD move(vector.z_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(z_ = z_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, z_ = z_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(pack.x);
		y_ = y_ - _STD move(pack.y);
		z_ = z_ - _STD move(pack.z);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(z_ *= scl))
		-> decltype(z_ *= scl, _STD declval<LRefEucVector>()) {
		x_ *= scl;
		y_ *= scl;
		z_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(z_ = z_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, z_ = z_ * scl, _STD declval<LRefEucVector>()) {
		x_ = x_ * scl;
		y_ = y_ * scl;
		z_ = z_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(z_ /= scl))
		-> decltype(z_ /= scl, _STD declval<LRefEucVector>()) {
		x_ /= scl;
		y_ /= scl;
		z_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(z_ = z_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, z_ = z_ / scl, _STD declval<LRefEucVector>()) {
		x_ = x_ / scl;
		y_ = y_ / scl;
		z_ = z_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return  y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType z() noexcept { return z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType z() const noexcept { return z_; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xx() const noexcept { return { x_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xy() const noexcept { return { x_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xz() const noexcept { return { x_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yx() const noexcept { return { y_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yy() const noexcept { return { y_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yz() const noexcept { return { y_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zx() const noexcept { return { z_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zy() const noexcept { return { z_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zz() const noexcept { return { z_,z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxx() const noexcept { return { x_,x_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxy() const noexcept { return { x_,x_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxz() const noexcept { return { x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyx() const noexcept { return { x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyy() const noexcept { return { x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyz() const noexcept { return { x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzx() const noexcept { return { x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzy() const noexcept { return { x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzz() const noexcept { return { x_, z_, z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxx() const noexcept { return { y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxy() const noexcept { return { y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxz() const noexcept { return { y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyx() const noexcept { return { y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyy() const noexcept { return { y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyz() const noexcept { return { y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzx() const noexcept { return { y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzy() const noexcept { return { y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzz() const noexcept { return { y_, z_, z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxx() const noexcept { return { z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxy() const noexcept { return { z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxz() const noexcept { return { z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyx() const noexcept { return { z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyy() const noexcept { return { z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyz() const noexcept { return { z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzx() const noexcept { return { z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzy() const noexcept { return { z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzz() const noexcept { return { z_, z_, z_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { x_, y_, z_ }; }
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		x_ = 0;
		y_ = 0;
		z_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y, class Z>
	EUCVECTORINLINE auto set(X&& x, Y&& y, Z&& z) noexcept(noexcept(x_ = _STD forward<X>(x)) && noexcept(y_ = _STD forward<Y>(y)) && noexcept(z_ = _STD forward<Z>(z)))
		-> meta::is_type_t<decltype(x_ = _STD forward<X>(x), y_ = _STD forward<Y>(y), z_ = _STD forward<Z>(z)), void> {
		x_ = _STD forward<X>(x);
		y_ = _STD forward<Y>(y);
		z_ = _STD forward<Z>(z);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanCmplVector3<T>& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* vector.x_) + ElemType(y_ * vector.y_))))
		-> decltype(ElemType(ElemType(x_* vector.x_) + ElemType(y_ * vector.y_)), x_ + x_) {
		return ElemType(x_ * vector.x_) + ElemType(y_ * vector.y_) + ElemType(z_ * vector.z_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanCmplVector3<T>&& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* _STD move(vector.x_)) + y_ * ElemType(_STD move(vector.y_)))))
		-> decltype(ElemType(ElemType(x_* _STD move(vector.x_)) + ElemType(y_ * _STD move(vector.y_))), x_ + x_) {
		return ElemType(x_ * _STD move(vector.x_)) + ElemType(y_ * _STD move(vector.y_)) + ElemType(z_ * _STD move(vector.z_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y)))))
		-> decltype(ElemType(ElemType(x_* _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y))), x_ + x_) {
		return ElemType(x_ * _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y)) + ElemType(z_ * _STD move(pack.z));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return ElemType(x_ * vector.x_) + ElemType(y_ * vector.y_) + ElemType(z_ * vector.z_);
	}
	/*
		@brief

			Calculate the norm.

			out = root((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt(ElemType(x_ * vector.x_) + ElemType(y_ * vector.y_) + ElemType(z_ * vector.z_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{z_ / eucnorm<T>()}))
		->decltype(Packer<T>{z_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { x_ / norm, y_ / norm, z_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}
	/*
		@brief

			Calculate the cross product.

	*/
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(const EuclideanCmplVector3<T>& right) const noexcept(
			noexcept(Packer<decltype(z_* right.z_ - z_ * right.z_)>{ z_* right.z_ - z_ * right.z_ }))
		-> decltype(Packer<decltype(z_* right.z_ - z_ * right.z_)>{ z_* right.z_ - z_ * right.z_ }) {
		return { y_ * right.z_ - z_ * right.y_, z_ * right.x_ - x_ * right.z_, x_ * right.y_ - y_ * right.x_ };
	}
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(EuclideanCmplVector3<T>&& right) const noexcept(
			noexcept(Packer<decltype(z_* _STD move(right.z_) - z_ * _STD move(right.z_))>{ z_* _STD move(right.z_) - z_ * _STD move(right.z_) }))
		-> decltype(Packer<decltype(z_* _STD move(right.z_) - z_ * _STD move(right.z_))>{ z_* _STD move(right.z_) - z_ * _STD move(right.z_) }) {
		return { y_ * _STD move(right.z_) - z_ * _STD move(right.y_), z_ * _STD move(right.x_) - x_ * _STD move(right.z_), x_ * _STD move(right.y_) - y_ * _STD move(right.x_) };
	}
	template<class T>
	EUCNODISCARD_MSG("The cross product calculation results are being ignored, which could suggest an unintended call.")
		EUCVECTORINLINE auto cross(RRefPacker<T> right) const noexcept(
			noexcept(Packer<decltype(z_* _STD move(right.z) - z_ * _STD move(right.z))>{ z_* _STD move(right.z) - z_ * _STD move(right.z) }))
		-> decltype(Packer<decltype(z_* _STD move(right.z) - z_ * _STD move(right.z))>{ z_* _STD move(right.z) - z_ * _STD move(right.z) }) {
		return { y_ * _STD move(right.z) - z_ * _STD move(right.y), z_ * _STD move(right.x) - x_ * _STD move(right.z), x_ * _STD move(right.y) - y_ * _STD move(right.x) };
	}
};

/*
	D4.
*/
template<class E = float>
struct EuclideanCmplVector4 final
	: private meta::evd_euc_vec {
protected:

	static constexpr size_t EucD = 4;
	static constexpr size_t RefSize = sizeof(E) * EucD;

	using ElemType = E;
	using EucVector = EuclideanCmplVector4<ElemType>;
	using ConstEucVector = _STD add_const_t<EucVector>;
	using LRefEucVector = _STD add_lvalue_reference_t<EucVector>;
	using LRefConstEucVector = _STD add_lvalue_reference_t<ConstEucVector>;
	using RRefEucVector = _STD add_rvalue_reference_t<EucVector>;
	using ConstElemType = _STD add_const_t<E>;
	using LRefElemType = _STD add_lvalue_reference_t<E>;
	using LRefConstElemType = _STD add_lvalue_reference_t<ConstElemType>;
	using RRefElemType = _STD add_rvalue_reference_t<E>;

	template<class R>
	using Packer = detail::ResultPacker_4<R>;
	template<class R>
	using SubPacker2 = detail::ResultPacker_2<R>;
	template<class R>
	using SubPacker3 = detail::ResultPacker_3<R>;
	template<class R>
	using RRefPacker = Packer<R>&&;

	template<class FE>
	friend struct EuclideanCmplVector4;

	static_assert(!_STD is_reference_v		<ElemType>, "Reference types arent allowed");
	static_assert(!_STD is_const_v			<ElemType>, "Member type must be mutable");
	static_assert(_STD is_constructible_v	<ElemType>, "Input type must be constructible without arguments");

public:

	ElemType x_, y_, z_, w_;

	/*
		Constructors.
	*/
	EuclideanCmplVector4() noexcept(_STD is_nothrow_constructible_v<ElemType>)
		: x_()
		, y_()
		, z_()
		, w_()
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ConstElemType>> = 0>
	EuclideanCmplVector4(const EuclideanCmplVector4& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ConstElemType>)
		: x_(vector.x_)
		, y_(vector.y_)
		, z_(vector.z_)
		, w_(vector.w_)
	{}

	template<class T = ElemType, meta::co_inst_if_t<T, ElemType, _STD is_constructible_v<ElemType, ElemType>> = 0>
	EuclideanCmplVector4(EuclideanCmplVector4&& vector) noexcept(_STD is_nothrow_constructible_v<ElemType, ElemType>)
		: x_(_STD move(vector.x_))
		, y_(_STD move(vector.y_))
		, z_(_STD move(vector.z_))
		, w_(_STD move(vector.w_))
	{}

	template<class T, meta::if_t<_STD is_constructible_v<ElemType, T>> = 0>
	EuclideanCmplVector4(RRefPacker<T> pack) noexcept(_STD is_nothrow_constructible_v<ElemType, T>)
		: x_(_STD move(pack.x))
		, y_(_STD move(pack.y))
		, z_(_STD move(pack.z))
		, w_(_STD move(pack.w))
	{}

	template<class X, class Y, class Z, class W, meta::if_t<meta::is_constructible_anynum_param_v<ElemType, X, Y, Z, W>> = 0>
	EuclideanCmplVector4(X&& x, Y&& y, Z&& z, W&& w) noexcept(_STD is_nothrow_constructible_v<ElemType, X> && _STD is_nothrow_constructible_v<ElemType, Y> && _STD is_nothrow_constructible_v<ElemType, Z> && _STD is_nothrow_constructible_v<ElemType, W>)
		: x_(_STD forward<X>(x))
		, y_(_STD forward<Y>(y))
		, z_(_STD forward<Z>(z))
		, w_(_STD forward<W>(w))
	{}

	/*
		Binary Operators.
	*/
	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(const EuclideanCmplVector4<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + vector.w_)>>{w_ + vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(w_ + vector.w_)>>{w_ + vector.w_}) {
		return { x_ + vector.x_, y_ + vector.y_, z_ + vector.z_, w_ + vector.w_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(EuclideanCmplVector4<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + _STD move(vector.w_))>>{w_ + _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(w_ + _STD move(vector.w_))>>{w_ + _STD move(vector.w_)}) {
		return { x_ + _STD move(vector.x_), y_ + _STD move(vector.y_), z_ + _STD move(vector.z_), w_ + _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE auto operator+(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ + _STD move(pack.x))>>{w_ + _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(w_ + _STD move(pack.x))>>{w_ + _STD move(pack.x)}) {
		return { x_ + _STD move(pack.x), y_ + _STD move(pack.y), z_ + _STD move(pack.z), w_ + _STD move(pack.w) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.w_)>>{_STD move(pack.x) + vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + vector.w_)>>{_STD move(pack.x) + vector.w_}) {
		return { _STD move(pack.x) + vector.x_ , _STD move(pack.y) + vector.y_ , _STD move(pack.z) + vector.z_ , _STD move(pack.w) + vector.w_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the addition is being ignored. If you intend to modify the lvalue, please use [+=] instead.")
		EUCVECTORINLINE static auto operator+(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.w_))>>{_STD move(pack.x) + _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) + _STD move(vector.w_))>>{_STD move(pack.x) + _STD move(vector.w_)}) {
		return { _STD move(pack.x) + _STD move(vector.x_), _STD move(pack.y) + _STD move(vector.y_),  _STD move(pack.z) + _STD move(vector.z_), _STD move(pack.w) + _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(const EuclideanCmplVector4<T>& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - vector.w_)>>{w_ - vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(w_ - vector.w_)>>{w_ - vector.w_}) {
		return { x_ - vector.x_, y_ - vector.y_, z_ - vector.z_, w_ - vector.w_ };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(EuclideanCmplVector4<T>&& vector) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - _STD move(vector.w_))>>{w_ - _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(w_ - _STD move(vector.w_))>>{w_ - _STD move(vector.w_)}) {
		return { x_ - _STD move(vector.x_), y_ - _STD move(vector.y_), z_ - _STD move(vector.z_), w_ - _STD move(vector.w_) };
	}

	template<class T>
	EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE auto operator-(RRefPacker<T> pack) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ - _STD move(pack.x))>>{w_ - _STD move(pack.x)}))
		->decltype(Packer<meta::no_ref<decltype(w_ - _STD move(pack.x))>>{w_ - _STD move(pack.x)}) {
		return { x_ - _STD move(pack.x), y_ - _STD move(pack.y), z_ - _STD move(pack.z), w_ - _STD move(pack.w) };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, LRefConstEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.w_)>>{_STD move(pack.x) - vector.w_}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - vector.w_)>>{_STD move(pack.x) - vector.w_}) {
		return { _STD move(pack.x) - vector.x_ , _STD move(pack.y) - vector.y_ , _STD move(pack.z) - vector.z_ , _STD move(pack.w) - vector.w_ };
	}

	template<class T>
	friend EUCNODISCARD_MSG("The result of the subtraction is being ignored. If you intend to modify the lvalue, please use [-=] instead.")
		EUCVECTORINLINE static auto operator-(RRefPacker<T> pack, RRefEucVector vector) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.w_))>>{_STD move(pack.x) - _STD move(vector.w_)}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(pack.x) - _STD move(vector.w_))>>{_STD move(pack.x) - _STD move(vector.w_)}) {
		return { _STD move(pack.x) - _STD move(vector.x_), _STD move(pack.y) - _STD move(vector.y_),  _STD move(pack.z) - _STD move(vector.z_), _STD move(pack.w) - _STD move(vector.w_) };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE auto operator*(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_* scl)>>{w_* scl}))
		->decltype(Packer<meta::no_ref<decltype(w_* scl)>>{w_* scl}) {
		return { x_ * scl, y_ * scl, z_ * scl, w_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, LRefConstEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(right.w_* scl)>>{right.w_* scl}))
		->decltype(Packer<meta::no_ref<decltype(right.w_* scl)>>{right.w_* scl}) {
		return { right.x_ * scl,  right.y_ * scl, right.z_ * scl, right.w_ * scl };
	}

	template<class S>
	friend EUCNODISCARD_MSG("The result of the multiplication is being ignored.If you intend to modify the lvalue, please use [*=] instead.")
		EUCVECTORINLINE static auto operator*(S&& scl, RRefEucVector right) noexcept(noexcept(Packer<meta::no_ref<decltype(_STD move(right.w_)* scl)>>{_STD move(right.w_)* scl}))
		->decltype(Packer<meta::no_ref<decltype(_STD move(right.w_)* scl)>>{_STD move(right.w_)* scl}) {
		return { _STD move(right.x_) * scl, _STD move(right.y_) * scl, _STD move(right.z_) * scl, _STD move(right.w_) * scl };
	}

	template<class S>
	EUCNODISCARD_MSG("The result of the division is being ignored.If you intend to modify the lvalue, please use [/=] instead.")
		EUCVECTORINLINE auto operator/(S&& scl) const noexcept(noexcept(Packer<meta::no_ref<decltype(w_ / scl)>>{w_ / scl}))
		->decltype(Packer<meta::no_ref<decltype(w_ / scl)>>{w_ / scl}) {
		return { x_ / scl, y_ / scl,  z_ / scl, w_ / scl };
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(const EuclideanCmplVector4<T>& right) const noexcept(noexcept(bool(w_ == right.w_)))
		-> meta::no_ref<decltype(bool(w_ == right.w_))> {
		return (x_ == right.x_) && (y_ == right.y_) && (z_ == right.z_) && (w_ == right.w_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(EuclideanCmplVector4<T>&& right) const noexcept(noexcept(bool(w_ == _STD move(right.w_))))
		-> meta::no_ref<decltype(bool(w_ == _STD move(right.w_)))> {
		return (x_ == _STD move(right.x_)) && (y_ == _STD move(right.y_)) && (z_ == _STD move(right.z_)) && (w_ == _STD move(right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator==(RRefPacker<T> right) const noexcept(noexcept(bool(w_ == _STD move(right.x))))
		-> meta::no_ref<decltype(bool(w_ == _STD move(right.x)))> {
		return (x_ == _STD move(right.x)) && (y_ == _STD move(right.y)) && (z_ == _STD move(right.z)) && (w_ == _STD move(right.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.w_ == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.w_ == _STD move(left.x)))> {
		return (right.x_ == _STD move(left.x)) && (right.y_ == _STD move(left.y)) && (right.z_ == _STD move(left.z)) && (right.w_ == _STD move(left.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator==(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.w_) == _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.w_) == _STD move(left.x)))> {
		return (_STD move(right.x_) == _STD move(left.x)) && (_STD move(right.y_) == _STD move(left.y)) && (_STD move(right.z_) == _STD move(left.z)) && (_STD move(right.w_) == _STD move(left.w));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector4<T>& right) const noexcept(noexcept(bool(w_ != right.w_)))
		-> meta::no_ref<decltype(bool(w_ != right.w_))> {
		return (x_ != right.x_) || (y_ != right.y_) || (z_ != right.z_) || (w_ != right.w_);
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector4<T>&& right) const noexcept(noexcept(bool(w_ != _STD move(right.w_))))
		-> meta::no_ref<decltype(bool(w_ != _STD move(right.w_)))> {
		return (x_ != _STD move(right.x_)) || (y_ != _STD move(right.y_)) || (z_ != _STD move(right.z_)) || (w_ != _STD move(right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(w_ != _STD move(right.x))))
		-> meta::no_ref<decltype(bool(w_ != _STD move(right.x)), _STD declval<bool>())> {
		return (x_ != _STD move(right.x)) || (y_ != _STD move(right.y)) || (z_ != _STD move(right.z)) || (w_ != _STD move(right.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(right.w_ != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(right.w_ != _STD move(left.x)), _STD declval<bool>())> {
		return (right.x_ != _STD move(left.x)) || (right.y_ != _STD move(left.y)) || (right.z_ != _STD move(left.z)) || (right.w_ != _STD move(left.w));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(_STD move(right.w_) != _STD move(left.x))))
		-> meta::no_ref<decltype(bool(_STD move(right.w_) != _STD move(left.x)), _STD declval<bool>())> {
		return (_STD move(right.x_) != _STD move(left.x)) || (_STD move(right.y_) != _STD move(left.y)) || (_STD move(right.z_) != _STD move(left.z)) || (_STD move(right.w_) != _STD move(left.w));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(const EuclideanCmplVector4<T>& right) const noexcept(noexcept(bool(!(w_ == right.w_))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(right)>>(), bool(!(w_ == right.w_)), _STD declval<bool>())> {
		return (!(x_ == right.x_)) || (!(y_ == right.y_)) || (!(z_ == right.z_)) || (!(w_ == right.w_));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(EuclideanCmplVector4<T>&& right) const noexcept(noexcept(bool(!(w_ == _STD move(right.w_)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(w_ == _STD move(right.w_))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x_))) || (!(y_ == _STD move(right.y_))) || (!(z_ == _STD move(right.z_))) || (!(w_ == _STD move(right.w_)));
	}

	template<class T>
	EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE auto operator!=(RRefPacker<T> right) const noexcept(noexcept(bool(!(w_ == _STD move(right.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(*this), decltype(_STD move(right))>>(), bool(!(w_ == _STD move(right.x))), _STD declval<bool>())> {
		return (!(x_ == _STD move(right.x))) || (!(y_ == _STD move(right.y))) || (!(z_ == _STD move(right.z))) || (!(w_ == _STD move(right.w)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, LRefConstEucVector right) noexcept(noexcept(bool(!(right.w_ == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(right)>>(), bool(!(right.w_ == _STD move(left.x))), _STD declval<bool>())> {
		return (!(right.x_ == _STD move(left.x))) || (!(right.y_ == _STD move(left.y))) || (!(right.z_ == _STD move(left.z))) || (!(right.w_ == _STD move(left.w)));
	}

	template<class T>
	friend EUCNODISCARD_MSG("The return value of the comparison operator is being ignored")
		EUCVECTORINLINE static auto operator!=(RRefPacker<T> left, RRefEucVector right) noexcept(noexcept(bool(!(_STD move(right.x) == _STD move(left.x)))))
		-> meta::no_ref<decltype(meta::cnd_type<!meta::is_not_equal_v<decltype(_STD move(left)), decltype(_STD move(right))>>(), bool(!(_STD move(right.w_)) == _STD move(left.x)), _STD declval<bool>())> {
		return (!(_STD move(right.x_) == _STD move(left.x))) || (!(_STD move(right.y_) == _STD move(left.y))) || (!(_STD move(right.z_) == _STD move(left.z))) || (!(_STD move(right.w_) == _STD move(left.w)));
	}
	/*
		Unary Operators.
	*/
	EUCNODISCARD EUCVECTORINLINE EuclideanCmplVector4<ElemType>& operator+() noexcept {
		return *this;
	}
	EUCNODISCARD EUCVECTORINLINE const EuclideanCmplVector4<ElemType>& operator+() const noexcept {
		return *this;
	}

	template<class T = ElemType>
	EUCNODISCARD EUCVECTORINLINE auto operator-() const noexcept(noexcept(EuclideanCmplVector4<T>() * -1))
		->decltype(EuclideanCmplVector4<T>() * -1) {
		return { x_ * -1, y_ * -1,  z_ * -1, w_ * -1 };
	}
	/*
		Assignment Operators.
	*/
	template<class T>
	EUCVECTORINLINE auto operator=(const EuclideanCmplVector4<T>& vector) & noexcept(noexcept(w_ = vector.w_))
		-> decltype(w_ = vector.w_, _STD declval<LRefEucVector>()) {
		x_ = vector.x_;
		y_ = vector.y_;
		z_ = vector.z_;
		w_ = vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(EuclideanCmplVector4<T>&& vector) & noexcept(noexcept(w_ = _STD move(vector.w_)))
		-> decltype(w_ = _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		x_ = _STD move(vector.x_);
		y_ = _STD move(vector.y_);
		z_ = _STD move(vector.z_);
		w_ = _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator=(RRefPacker<T> pack) & noexcept(noexcept(w_ = _STD move(pack.x)))
		-> decltype(w_ = _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = _STD move(pack.x);
		y_ = _STD move(pack.y);
		z_ = _STD move(pack.z);
		w_ = _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector4<T>& vector) & noexcept(noexcept(w_ += vector.w_))
		-> decltype(w_ += vector.w_, _STD declval<LRefEucVector>()) {
		x_ += vector.x_;
		y_ += vector.y_;
		z_ += vector.z_;
		w_ += vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector4<T>&& vector) & noexcept(noexcept(w_ += _STD move(vector.w_)))
		-> decltype(w_ += _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		x_ += _STD move(vector.x_);
		y_ += _STD move(vector.y_);
		z_ += _STD move(vector.z_);
		w_ += _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(w_ += _STD move(pack.x)))
		-> decltype(w_ += _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ += _STD move(pack.x);
		y_ += _STD move(pack.y);
		z_ += _STD move(pack.z);
		w_ += _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(const EuclideanCmplVector4<T>& vector) & noexcept(noexcept(w_ = w_ + vector.w_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, w_ = w_ + vector.w_, _STD declval<LRefEucVector>()) {
		x_ = x_ + vector.x_;
		y_ = y_ + vector.y_;
		z_ = z_ + vector.z_;
		w_ = w_ + vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(EuclideanCmplVector4<T>&& vector) & noexcept(noexcept(w_ = w_ + _STD move(vector.w_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, w_ = w_ + _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(vector.x_);
		y_ = y_ + _STD move(vector.y_);
		z_ = z_ + _STD move(vector.z_);
		w_ = w_ + _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator+=(RRefPacker<T> pack) & noexcept(noexcept(w_ = w_ + _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, w_ = w_ + _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ + _STD move(pack.x);
		y_ = y_ + _STD move(pack.y);
		z_ = z_ + _STD move(pack.z);
		w_ = w_ + _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector4<T>& vector) & noexcept(noexcept(w_ -= vector.w_))
		-> decltype(w_ -= vector.w_, _STD declval<LRefEucVector>()) {
		x_ -= vector.x_;
		y_ -= vector.y_;
		z_ -= vector.z_;
		w_ -= vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector4<T>&& vector) & noexcept(noexcept(w_ -= _STD move(vector.w_)))
		-> decltype(w_ -= _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(vector.x_);
		y_ -= _STD move(vector.y_);
		z_ -= _STD move(vector.z_);
		w_ -= _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(w_ -= _STD move(pack.x)))
		-> decltype(w_ -= _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ -= _STD move(pack.x);
		y_ -= _STD move(pack.y);
		z_ -= _STD move(pack.z);
		w_ -= _STD move(pack.w);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(const EuclideanCmplVector4<T>& vector) & noexcept(noexcept(w_ = w_ - vector.w_))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(vector)>>, w_ = w_ - vector.w_, _STD declval<LRefEucVector>()) {
		x_ = x_ - vector.x_;
		y_ = y_ - vector.y_;
		z_ = z_ - vector.z_;
		w_ = w_ - vector.w_;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(EuclideanCmplVector4<T>&& vector) & noexcept(noexcept(w_ = w_ - _STD move(vector.w_)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(vector))>>, w_ = w_ - _STD move(vector.w_), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(vector.x_);
		y_ = y_ - _STD move(vector.y_);
		z_ = z_ - _STD move(vector.z_);
		w_ = w_ - _STD move(vector.w_);
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator-=(RRefPacker<T> pack) & noexcept(noexcept(w_ = w_ - _STD move(pack.x)))
		-> decltype(meta::cnd_type<!meta::is_add_equal_v<decltype(*this), decltype(_STD move(pack))>>, w_ = w_ - _STD move(pack.x), _STD declval<LRefEucVector>()) {
		x_ = x_ - _STD move(pack.x);
		y_ = y_ - _STD move(pack.y);
		z_ = z_ - _STD move(pack.z);
		w_ = w_ - _STD move(pack.w);
		return *this;
	}


	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(w_ *= scl))
		-> decltype(w_ *= scl, _STD declval<LRefEucVector>()) {
		x_ *= scl;
		y_ *= scl;
		z_ *= scl;
		w_ *= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator*=(T&& scl) & noexcept(noexcept(w_ = w_ * scl))
		-> decltype(meta::cnd_type<!meta::is_mul_equal_v<decltype(*this), decltype(scl)>>, w_ = w_ * scl, _STD declval<LRefEucVector>()) {
		x_ = x_ * scl;
		y_ = y_ * scl;
		z_ = z_ * scl;
		w_ = w_ * scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(w_ /= scl))
		-> decltype(w_ /= scl, _STD declval<LRefEucVector>()) {
		x_ /= scl;
		y_ /= scl;
		z_ /= scl;
		w_ /= scl;
		return *this;
	}

	template<class T>
	EUCVECTORINLINE auto operator/=(T&& scl) & noexcept(noexcept(w_ = w_ / scl))
		-> decltype(meta::cnd_type<!meta::is_div_equal_v<decltype(*this), decltype(scl)>>, w_ = w_ / scl, _STD declval<LRefEucVector>()) {
		x_ = x_ / scl;
		y_ = y_ / scl;
		z_ = z_ / scl;
		w_ = w_ / scl;
		return *this;
	}

	/*
		Client Function.
	*/
	/*
		@brief

			Get dimension.

	*/
	EUCNODISCARD_MSG("The acquisition of dimensionality is disregarded. It is possible that this is an unintended call.")
		EUCVECTORINLINE constexpr size_t dimension() { return EucD; }
	/*
		@brief

			Obtaining a vector of elements or lower dimensions.

	*/
#ifdef _MSC_VER
#pragma region 
#endif
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType x() noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType x() const noexcept { return x_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType y() noexcept { return y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType y() const noexcept { return  y_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType z() noexcept { return z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType z() const noexcept { return z_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefElemType w() noexcept { return w_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE constexpr LRefConstElemType w() const noexcept { return w_; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xx() const noexcept { return { x_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xy() const noexcept { return { x_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xz() const noexcept { return { x_,z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> xw() const noexcept { return { x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yx() const noexcept { return { y_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yy() const noexcept { return { y_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yz() const noexcept { return { y_,z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> yw() const noexcept { return { y_,w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zx() const noexcept { return { z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zy() const noexcept { return { z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zz() const noexcept { return { z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> zw() const noexcept { return { z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wx() const noexcept { return { w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wy() const noexcept { return { w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> wz() const noexcept { return { w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker2<ElemType> ww() const noexcept { return { w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxx() const noexcept { return { x_,x_,x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxy() const noexcept { return { x_,x_,y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxz() const noexcept { return { x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xxw() const noexcept { return { x_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyx() const noexcept { return { x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyy() const noexcept { return { x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyz() const noexcept { return { x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xyw() const noexcept { return { x_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzx() const noexcept { return { x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzy() const noexcept { return { x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzz() const noexcept { return { x_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xzw() const noexcept { return { x_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwx() const noexcept { return { x_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwy() const noexcept { return { x_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xwz() const noexcept { return { x_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> xww() const noexcept { return { x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxx() const noexcept { return { y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxy() const noexcept { return { y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxz() const noexcept { return { y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yxw() const noexcept { return { y_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyx() const noexcept { return { y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyy() const noexcept { return { y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyz() const noexcept { return { y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yyw() const noexcept { return { y_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzx() const noexcept { return { y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzy() const noexcept { return { y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzz() const noexcept { return { y_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yzw() const noexcept { return { y_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywx() const noexcept { return { y_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywy() const noexcept { return { y_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> ywz() const noexcept { return { y_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> yww() const noexcept { return { y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxx() const noexcept { return { z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxy() const noexcept { return { z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxz() const noexcept { return { z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zxw() const noexcept { return { z_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyx() const noexcept { return { z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyy() const noexcept { return { z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyz() const noexcept { return { z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zyw() const noexcept { return { z_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzx() const noexcept { return { z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzy() const noexcept { return { z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzz() const noexcept { return { z_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zzw() const noexcept { return { z_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwx() const noexcept { return { z_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwy() const noexcept { return { z_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zwz() const noexcept { return { z_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> zww() const noexcept { return { z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxx() const noexcept { return { w_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxy() const noexcept { return { w_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxz() const noexcept { return { w_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wxw() const noexcept { return { w_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyx() const noexcept { return { w_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyy() const noexcept { return { w_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyz() const noexcept { return { w_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wyw() const noexcept { return { w_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzx() const noexcept { return { w_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzy() const noexcept { return { w_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzz() const noexcept { return { w_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wzw() const noexcept { return { w_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwx() const noexcept { return { w_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwy() const noexcept { return { w_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> wwz() const noexcept { return { w_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE SubPacker3<ElemType> www() const noexcept { return { w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxx() const noexcept { return { x_, x_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxy() const noexcept { return { x_, x_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxz() const noexcept { return { x_, x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxxw() const noexcept { return { x_, x_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyx() const noexcept { return { x_, x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyy() const noexcept { return { x_, x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyz() const noexcept { return { x_, x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxyw() const noexcept { return { x_, x_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzx() const noexcept { return { x_, x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzy() const noexcept { return { x_, x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzz() const noexcept { return { x_, x_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxzw() const noexcept { return { x_, x_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwx() const noexcept { return { x_, x_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwy() const noexcept { return { x_, x_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxwz() const noexcept { return { x_, x_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xxww() const noexcept { return { x_, x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxx() const noexcept { return { x_, y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxy() const noexcept { return { x_, y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxz() const noexcept { return { x_, y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyxw() const noexcept { return { x_, y_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyx() const noexcept { return { x_, y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyy() const noexcept { return { x_, y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyz() const noexcept { return { x_, y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyyw() const noexcept { return { x_, y_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzx() const noexcept { return { x_, y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzy() const noexcept { return { x_, y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyzz() const noexcept { return { x_, y_, z_, z_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywx() const noexcept { return { x_, y_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywy() const noexcept { return { x_, y_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xywz() const noexcept { return { x_, y_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xyww() const noexcept { return { x_, y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxx() const noexcept { return { x_, z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxy() const noexcept { return { x_, z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxz() const noexcept { return { x_, z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzxw() const noexcept { return { x_, z_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyx() const noexcept { return { x_, z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyy() const noexcept { return { x_, z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyz() const noexcept { return { x_, z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzyw() const noexcept { return { x_, z_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzx() const noexcept { return { x_, z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzy() const noexcept { return { x_, z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzz() const noexcept { return { x_, z_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzzw() const noexcept { return { x_, z_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwx() const noexcept { return { x_, z_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwy() const noexcept { return { x_, z_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzwz() const noexcept { return { x_, z_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xzww() const noexcept { return { x_, z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxx() const noexcept { return { x_, w_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxy() const noexcept { return { x_, w_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxz() const noexcept { return { x_, w_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwxw() const noexcept { return { x_, w_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyx() const noexcept { return { x_, w_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyy() const noexcept { return { x_, w_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyz() const noexcept { return { x_, w_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwyw() const noexcept { return { x_, w_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzx() const noexcept { return { x_, w_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzy() const noexcept { return { x_, w_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzz() const noexcept { return { x_, w_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwzw() const noexcept { return { x_, w_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwx() const noexcept { return { x_, w_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwy() const noexcept { return { x_, w_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwwz() const noexcept { return { x_, w_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> xwww() const noexcept { return { x_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxx() const noexcept { return { y_, x_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxy() const noexcept { return { y_, x_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxz() const noexcept { return { y_, x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxxw() const noexcept { return { y_, x_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyx() const noexcept { return { y_, x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyy() const noexcept { return { y_, x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyz() const noexcept { return { y_, x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxyw() const noexcept { return { y_, x_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzx() const noexcept { return { y_, x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzy() const noexcept { return { y_, x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzz() const noexcept { return { y_, x_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxzw() const noexcept { return { y_, x_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwx() const noexcept { return { y_, x_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwy() const noexcept { return { y_, x_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxwz() const noexcept { return { y_, x_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yxww() const noexcept { return { y_, x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxx() const noexcept { return { y_, y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxy() const noexcept { return { y_, y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxz() const noexcept { return { y_, y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyxw() const noexcept { return { y_, y_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyx() const noexcept { return { y_, y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyy() const noexcept { return { y_, y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyz() const noexcept { return { y_, y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyyw() const noexcept { return { y_, y_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzx() const noexcept { return { y_, y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzy() const noexcept { return { y_, y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzz() const noexcept { return { y_, y_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyzw() const noexcept { return { y_, y_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywx() const noexcept { return { y_, y_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywy() const noexcept { return { y_, y_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yywz() const noexcept { return { y_, y_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yyww() const noexcept { return { y_, y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxx() const noexcept { return { y_, z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxy() const noexcept { return { y_, z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxz() const noexcept { return { y_, z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzxw() const noexcept { return { y_, z_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyx() const noexcept { return { y_, z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyy() const noexcept { return { y_, z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyz() const noexcept { return { y_, z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzyw() const noexcept { return { y_, z_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzx() const noexcept { return { y_, z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzy() const noexcept { return { y_, z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzz() const noexcept { return { y_, z_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzzw() const noexcept { return { y_, z_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwx() const noexcept { return { y_, z_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwy() const noexcept { return { y_, z_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzwz() const noexcept { return { y_, z_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> yzww() const noexcept { return { y_, z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxx() const noexcept { return { y_, w_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxy() const noexcept { return { y_, w_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxz() const noexcept { return { y_, w_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywxw() const noexcept { return { y_, w_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyx() const noexcept { return { y_, w_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyy() const noexcept { return { y_, w_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyz() const noexcept { return { y_, w_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywyw() const noexcept { return { y_, w_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzx() const noexcept { return { y_, w_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzy() const noexcept { return { y_, w_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzz() const noexcept { return { y_, w_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywzw() const noexcept { return { y_, w_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwx() const noexcept { return { y_, w_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwy() const noexcept { return { y_, w_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywwz() const noexcept { return { y_, w_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> ywww() const noexcept { return { y_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxx() const noexcept { return { z_, x_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxy() const noexcept { return { z_, x_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxz() const noexcept { return { z_, x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxxw() const noexcept { return { z_, x_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyx() const noexcept { return { z_, x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyy() const noexcept { return { z_, x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyz() const noexcept { return { z_, x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxyw() const noexcept { return { z_, x_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzx() const noexcept { return { z_, x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzy() const noexcept { return { z_, x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzz() const noexcept { return { z_, x_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxzw() const noexcept { return { z_, x_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwx() const noexcept { return { z_, x_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwy() const noexcept { return { z_, x_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxwz() const noexcept { return { z_, x_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zxww() const noexcept { return { z_, x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxx() const noexcept { return { z_, y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxy() const noexcept { return { z_, y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxz() const noexcept { return { z_, y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyxw() const noexcept { return { z_, y_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyx() const noexcept { return { z_, y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyy() const noexcept { return { z_, y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyz() const noexcept { return { z_, y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyyw() const noexcept { return { z_, y_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzx() const noexcept { return { z_, y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzy() const noexcept { return { z_, y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzz() const noexcept { return { z_, y_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyzw() const noexcept { return { z_, y_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywx() const noexcept { return { z_, y_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywy() const noexcept { return { z_, y_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zywz() const noexcept { return { z_, y_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zyww() const noexcept { return { z_, y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxx() const noexcept { return { z_, z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxy() const noexcept { return { z_, z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxz() const noexcept { return { z_, z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzxw() const noexcept { return { z_, z_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyx() const noexcept { return { z_, z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyy() const noexcept { return { z_, z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyz() const noexcept { return { z_, z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzyw() const noexcept { return { z_, z_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzx() const noexcept { return { z_, z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzy() const noexcept { return { z_, z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzz() const noexcept { return { z_, z_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzzw() const noexcept { return { z_, z_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwx() const noexcept { return { z_, z_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwy() const noexcept { return { z_, z_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzwz() const noexcept { return { z_, z_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zzww() const noexcept { return { z_, z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxx() const noexcept { return { z_, w_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxy() const noexcept { return { z_, w_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxz() const noexcept { return { z_, w_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwxw() const noexcept { return { z_, w_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyx() const noexcept { return { z_, w_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyy() const noexcept { return { z_, w_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyz() const noexcept { return { z_, w_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwyw() const noexcept { return { z_, w_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzx() const noexcept { return { z_, w_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzy() const noexcept { return { z_, w_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzz() const noexcept { return { z_, w_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwzw() const noexcept { return { z_, w_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwx() const noexcept { return { z_, w_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwy() const noexcept { return { z_, w_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwwz() const noexcept { return { z_, w_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> zwww() const noexcept { return { z_, w_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxx() const noexcept { return { w_, x_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxy() const noexcept { return { w_, x_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxz() const noexcept { return { w_, x_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxxw() const noexcept { return { w_, x_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyx() const noexcept { return { w_, x_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyy() const noexcept { return { w_, x_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyz() const noexcept { return { w_, x_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxyw() const noexcept { return { w_, x_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzx() const noexcept { return { w_, x_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzy() const noexcept { return { w_, x_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzz() const noexcept { return { w_, x_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxzw() const noexcept { return { w_, x_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwx() const noexcept { return { w_, x_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwy() const noexcept { return { w_, x_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxwz() const noexcept { return { w_, x_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wxww() const noexcept { return { w_, x_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxx() const noexcept { return { w_, y_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxy() const noexcept { return { w_, y_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxz() const noexcept { return { w_, y_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyxw() const noexcept { return { w_, y_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyx() const noexcept { return { w_, y_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyy() const noexcept { return { w_, y_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyz() const noexcept { return { w_, y_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyyw() const noexcept { return { w_, y_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzx() const noexcept { return { w_, y_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzy() const noexcept { return { w_, y_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzz() const noexcept { return { w_, y_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyzw() const noexcept { return { w_, y_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywx() const noexcept { return { w_, y_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywy() const noexcept { return { w_, y_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wywz() const noexcept { return { w_, y_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wyww() const noexcept { return { w_, y_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxx() const noexcept { return { w_, z_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxy() const noexcept { return { w_, z_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxz() const noexcept { return { w_, z_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzxw() const noexcept { return { w_, z_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyx() const noexcept { return { w_, z_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyy() const noexcept { return { w_, z_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyz() const noexcept { return { w_, z_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzyw() const noexcept { return { w_, z_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzx() const noexcept { return { w_, z_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzy() const noexcept { return { w_, z_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzz() const noexcept { return { w_, z_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzzw() const noexcept { return { w_, z_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwx() const noexcept { return { w_, z_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwy() const noexcept { return { w_, z_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzwz() const noexcept { return { w_, z_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wzww() const noexcept { return { w_, z_, w_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxx() const noexcept { return { w_, w_, x_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxy() const noexcept { return { w_, w_, x_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxz() const noexcept { return { w_, w_, x_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwxw() const noexcept { return { w_, w_, x_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyx() const noexcept { return { w_, w_, y_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyy() const noexcept { return { w_, w_, y_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyz() const noexcept { return { w_, w_, y_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwyw() const noexcept { return { w_, w_, y_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzx() const noexcept { return { w_, w_, z_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzy() const noexcept { return { w_, w_, z_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzz() const noexcept { return { w_, w_, z_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwzw() const noexcept { return { w_, w_, z_, w_ }; }

	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwx() const noexcept { return { w_, w_, w_, x_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwy() const noexcept { return { w_, w_, w_, y_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwwz() const noexcept { return { w_, w_, w_, z_ }; }
	EUCNODISCARD_MSG("The reference has been discarded. This may be an unintended call.")
		EUCVECTORINLINE Packer<ElemType> wwww() const noexcept { return { w_, w_, w_, w_ }; }
	EUCNODISCARD_MSG("The output has been discarded. There may have been an unintended call.")
		EUCVECTORINLINE Packer<ElemType> get_pack() const noexcept { return { x_,y_, z_, w_ }; }
#ifdef _MSC_VER
#pragma endregion
#endif
	/*
		@brief

			Make this vector into a zero vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto zero_self() noexcept(noexcept(_STD declval<_STD add_lvalue_reference_t<T>>() = 0))
		-> meta::is_type_t<decltype(_STD declval<_STD add_lvalue_reference_t<T>>() = 0), void> {
		x_ = 0;
		y_ = 0;
		z_ = 0;
		w_ = 0;
	}
	/*
		@brief

			Set values for elements.

			(x = in1,y = in2, z = in3, ..., n = inn).

	*/
	template<class X, class Y, class Z, class W>
	EUCVECTORINLINE auto set(X&& x, Y&& y, Z&& z, W&& w) noexcept(noexcept(x_ = _STD forward<X>(x)) && noexcept(y_ = _STD forward<Y>(y)) && noexcept(z_ = _STD forward<Z>(z)) && noexcept(w_ = _STD forward<W>(w)))
		-> meta::is_type_t<decltype(x_ = _STD forward<X>(x), y_ = _STD forward<Y>(y), z_ = _STD forward<Z>(z), z_ = _STD forward<W>(w)), void> {
		x_ = _STD forward<X>(x);
		y_ = _STD forward<Y>(y);
		z_ = _STD forward<Z>(z);
		w_ = _STD forward<W>(w);
	}
	/*
		@brief

			Calculates the dot product.

			out = e1*re1 + e2*re2 + ... + en*ren.

	*/
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(const EuclideanCmplVector4<T>& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* vector.x_) + ElemType(y_ * vector.y_))))
		-> decltype(ElemType(ElemType(x_* vector.x_) + ElemType(y_ * vector.y_)), x_ + x_) {
		return ElemType(x_ * vector.x_) + ElemType(y_ * vector.y_) + ElemType(z_ * vector.z_) + ElemType(w_ * vector.w_);
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(EuclideanCmplVector4<T>&& vector) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* _STD move(vector.x_)) + y_ * ElemType(_STD move(vector.y_)))))
		-> decltype(ElemType(ElemType(x_* _STD move(vector.x_)) + ElemType(y_ * _STD move(vector.y_))), x_ + x_) {
		return ElemType(x_ * _STD move(vector.x_)) + ElemType(y_ * _STD move(vector.y_)) + ElemType(z_ * _STD move(vector.z_)) + ElemType(w_ * _STD move(vector.w_));
	}
	template<class T>
	EUCNODISCARD_MSG("The dot product calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto dot(RRefPacker<T> pack) const noexcept(noexcept(x_ + x_) && noexcept(ElemType(ElemType(x_* _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y)))))
		-> decltype(ElemType(ElemType(x_* _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y))), x_ + x_) {
		return ElemType(x_ * _STD move(pack.x)) + ElemType(y_ * _STD move(pack.y)) + ElemType(z_ * _STD move(pack.z)) + ElemType(w_ * _STD move(pack.w));
	}
	/*
		@brief

			Calculate the square of the norm.

			out = e1^2 + e2^2 + e3^2 +...+en^2.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm squared calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm_squared() const noexcept(noexcept(dot(_STD declval<LRefEucVector>())))
		-> decltype(dot(_STD declval<LRefEucVector>())) {
		return (x_ * x_) + (y_ * y_) + (z_ * z_) + (w_ * w_);
	}
	/*
		@brief

			Calculate the norm.

			out = ROOT((e1^2 + e2^2 + e3^2 +...+en^2)).

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The norm calculation results were ignored. This may be an unintended call.")
		EUCVECTORINLINE auto eucnorm() const noexcept(noexcept(_STD sqrt(eucnorm_squared<T>())))
		-> decltype(_STD sqrt(eucnorm_squared<T>())) {
		return _STD sqrt((x_ * x_) + (y_ * y_) + (z_ * z_) + (w_ * w_));
	}
	/*
		@brief

			Returns the result of normalizing this vector.

	*/
	template<class T = ElemType>
	EUCNODISCARD_MSG("The result of the normalization calculation was ignored. If you actually want to normalize this vector, use [normalize_self].")
		EUCVECTORINLINE auto normalize() const noexcept(noexcept(Packer<T>{w_ / eucnorm<T>()}))
		->decltype(Packer<T>{w_ / eucnorm<T>()}) {
		auto&& norm = eucnorm<T>();
		return { x_ / norm, y_ / norm, z_ / norm, w_ / norm };
	}
	/*
		@brief

			Normalize this vector.

	*/
	template<class T = ElemType>
	EUCVECTORINLINE auto normalize_self() noexcept(noexcept(_STD declval<LRefEucVector>() /= eucnorm<T>()))
		-> decltype(_STD declval<LRefEucVector>() /= eucnorm<T>()) {
		return *this /= eucnorm<T>();
	}

};

/*
	Vector Computed Result.
*/
template<class T>
using EucComputedResult1 = detail::ResultPacker_1<T>;
template<class T>
using EucComputedResult2 = detail::ResultPacker_2<T>;
template<class T>
using EucComputedResult3 = detail::ResultPacker_3<T>;
template<class T>
using EucComputedResult4 = detail::ResultPacker_4<T>;

using EucFloatComputedResult1 = detail::ResultPacker_1<float>;
using EucFloatComputedResult2 = detail::ResultPacker_2<float>;
using EucFloatComputedResult3 = detail::ResultPacker_3<float>;
using EucFloatComputedResult4 = detail::ResultPacker_4<float>;

using EucIntComputedResult1 = detail::ResultPacker_1<int>;
using EucIntComputedResult2 = detail::ResultPacker_2<int>;
using EucIntComputedResult3 = detail::ResultPacker_3<int>;
using EucIntComputedResult4 = detail::ResultPacker_4<int>;

using EucDoubleComputedResult1 = detail::ResultPacker_1<double>;
using EucDoubleComputedResult2 = detail::ResultPacker_2<double>;
using EucDoubleComputedResult3 = detail::ResultPacker_3<double>;
using EucDoubleComputedResult4 = detail::ResultPacker_4<double>;

/*
	Basic Vector type.
*/
using EucRecFloatVector1 = EuclideanVector1<float>;
using EucRecFloatVector2 = EuclideanRecVector2<float>;
using EucRecFloatVector3 = EuclideanRecVector3<float>;
using EucRecFloatVector4 = EuclideanRecVector4<float>;

using EucRecIntVector1 = EuclideanVector1<int>;
using EucRecIntVector2 = EuclideanRecVector2<int>;
using EucRecIntVector3 = EuclideanRecVector3<int>;
using EucRecIntVector4 = EuclideanRecVector4<int>;

using EucRecDoubleVector1 = EuclideanVector1<double>;
using EucRecDoubleVector2 = EuclideanRecVector2<double>;
using EucRecDoubleVector3 = EuclideanRecVector3<double>;
using EucRecDoubleVector4 = EuclideanRecVector4<double>;

using EucCmplFloatVector1 = EuclideanVector1<float>;
using EucCmplFloatVector2 = EuclideanCmplVector2<float>;
using EucCmplFloatVector3 = EuclideanCmplVector3<float>;
using EucCmplFloatVector4 = EuclideanCmplVector4<float>;

using EucCmplIntVector1 = EuclideanVector1<int>;
using EucCmplIntVector2 = EuclideanCmplVector2<int>;
using EucCmplIntVector3 = EuclideanCmplVector3<int>;
using EucCmplIntVector4 = EuclideanCmplVector4<int>;

using EucCmplDoubleVector1 = EuclideanVector1<double>;
using EucCmplDoubleVector2 = EuclideanCmplVector2<double>;
using EucCmplDoubleVector3 = EuclideanCmplVector3<double>;
using EucCmplDoubleVector4 = EuclideanCmplVector4<double>;

//name space end.
};
};

#endif 
