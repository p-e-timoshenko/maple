# Package: VectorCalcululsTools
# Version: 1.0
# Author: Timoshenko, Pavel Evgenjevich
# License: GNU Lesser General Public License v3 (LGPLv3)
# Date: 2014-08-07
# 
# Usages:
# > restart:
# > with(VectorCalculus):
# > with(VectorCalculusTools):
# > SetCoordinates('cartesian'[x,y,z])
#
# 1. Define VectorField
# 1.1> VF(e);
#  Vector[column](3, [e[x](x, y, z), e[y](x, y, z), e[z](x, y, z)])
# 1.2> VF(e[x]);
#  Vector[column](3, [e[x](x, y, z), 0, 0])
# 1.3> VF(e(u,v));
#  Vector[column](3, [e[x](u, v), e[y](u, v), e[z](u, v)])
# 1.4> VF(e(u,__,v));
#  Vector[column](3, [e[x](u, x, y, z, v), e[y](u, x, y, z, v), e[z](u, x, y, z, v)])
# 1.5> VF(e, t);
#  Vector[column](3, [e[x](x, y, z, t), e[y](x, y, z, t), e[z](x, y, z, t)])
# 1.6> VF(e(u, __, v), t);
#  Vector[column](3, [e[x](u, x, y, z, v, t), e[y](u, x, y, z, v, t), e[z](u, x, y, z, v, t)])
# 1.7> VF(e, t, coords = 'spherical'[u, v, w]);
#  Vector[column](3, [e[u](u, v, w, t), e[v](u, v, w, t), e[w](u, v, w, t)])
#
# 2. SubsVF, ZeroVF
# 2.1> Laplacian(VF(e)) = ZeroVF();
#  Vector[column](3, [
#    diff(diff(e[x](x, y, z), x), x)+diff(diff(e[x](x, y, z), y), y)+diff(diff(e[x](x, y, z), z), z),
#    diff(diff(e[y](x, y, z), x), x)+diff(diff(e[y](x, y, z), y), y)+diff(diff(e[y](x, y, z), z), z), 
#    diff(diff(e[z](x, y, z), x), x)+diff(diff(e[z](x, y, z), y), y)+diff(diff(e[z](x, y, z), z), z)]) 
#  = Vector[column](3, [0, 0, 0])
# 2.2> (2.2): SubsVF(VF(e)=VF(e[x]), %);
#  Vector[column](3, [diff(diff(e[x](x, y, z), x), x)+diff(diff(e[x](x, y, z), y), y)+diff(diff(e[x](x, y, z), z), z), 0, 0]) 
#  = Vector[column](3, [0, 0, 0])
#
# 3. VFEquationToList
# 3.1> Laplacian(VF(e)) = ZeroVF();
#    diff(diff(e[x](x, y, z), x), x)+diff(diff(e[x](x, y, z), y), y)+diff(diff(e[x](x, y, z), z), z),
#    diff(diff(e[y](x, y, z), x), x)+diff(diff(e[y](x, y, z), y), y)+diff(diff(e[y](x, y, z), z), z), 
#    diff(diff(e[z](x, y, z), x), x)+diff(diff(e[z](x, y, z), y), y)+diff(diff(e[z](x, y, z), z), z)]) 
#  = Vector[column](3, [0, 0, 0])
# 3.2> (3.2): VFEquationToList(%);
#  [
#    diff(e[x](x, y, z), x, x)+diff(e[x](x, y, z), y, y)+diff(e[x](x, y, z), z, z) = 0,
#    diff(e[y](x, y, z), x, x)+diff(e[y](x, y, z), y, y)+diff(e[y](x, y, z), z, z) = 0,
#    diff(e[z](x, y, z), x, x)+diff(e[z](x, y, z), y, y)+diff(e[z](x, y, z), z, z) = 0
#  ]
#
# 4. FormatCoordinates
# 4.1> FormatCoordinates()
#  cartesian[x, y, z]
# 4.2> FormatCoordinates('spherical')
#  spherical[r, phi, theta]
#
restart;
kernelopts('ASSERT' = true, 'assertlevel' = 2):
interface('showassumed' = 0):

module VectorCalculusTools ()
	description "VectorCalculus tools package";
	option
		package,
		`Copyright (C) Pavel Evgenjevich Timoshenko, 2014`;
		
	export
		FormatCoordinates,
		ZeroVF,
		VF,
		VFEquationToList,
		SubsVF;

	local
		ModuleLoad,
		ModuleUnload,
		_ParseCoefficientParameters,
		_FormatCoefficientParameters,
		_ApplyTemplate,
		_ApplyTemplateCoordinate,
		_ApplyTemplateCoordinates;
	
	FormatCoordinates := proc(
		coords::Or(name, indexed(name)) := VectorCalculus:-GetCoordinates(),
	$)::indexed(name);
		if type(coords, 'indexed(name)') then
			return coords;
		elif ormap(verify, [
			'bipolar', 'cardioid', 'cassinian', 'elliptic', 'hyperbolic', 
			'invcassinian', 'logarithmic', 'logcosh', 'parabolic',
			'rose', 'tangent'], coords, 'symbol') then
			return coords['u', 'v'];
		elif ormap(verify, ['cartesian2d', 'cartesian_2d'], coords, 'symbol') then
			return 'cartesian'['x', 'y'];
		elif verify('polar', coords, 'symbol') then
			return coords['r', 'theta'];
		elif ormap(verify, [
			'cartesian', 'cartesian3d', 'cartesian_3d'], coords, 'symbol') then
			return 'cartesian'['x', 'y', 'z'];			
		elif verify('cylindrical', coords, 'symbol') then
			return coords['r', 'theta', 'z'];
		elif verify('spherical', coords, 'symbol') then
			return coords['r', 'phi', 'theta'];
		elif ormap(verify, 
			['bipolarcylindrical', 'bispherical', 'cardioidal',
			'cardioidcylindrical', 'casscylindrical', 'conical',
			'ellcylindrical', 'hypercylindrical', 'invcasscylindrical',
			'logcylindrical', 'logcoshcylindrical',	'oblatespheroidal',
			'paraboloidal', 'paracylindrical',	'prolatespheroidal',
			'rosecylindrical', 'sixsphere', 'tangentcylindrical',
			'tangentsphere', 'toroidal'], coords, 'symbol') then
			return coords['u', 'v', 'w'];			
		else
			error "invalid coords: %1", coords;
		end if;
	end proc;
	protect('FormatCoordinates');

	ZeroVF := proc(
		{coords :: Or(name, indexed(name)) := VectorCalculus:-GetCoordinates()},
	$)::Vector[column];
		local
			_coords :: indexed(name);
			
		_coords := FormatCoordinates(coords);
		
		return VectorCalculus:-VectorField(convert(
			[0 $ nops(_coords)], 'Vector'), _coords); 
	end proc;
	protect('ZeroVF');

	VFEquationToList:= proc(
		eq :: equation(Vector[column]),
	$) :: list(equation);
		local
			leq :: list;

		leq := [lhs(eq), rhs(eq)];
		if not verify(op(map(
			VectorCalculus:-GetCoordinates, leq))) then
			error "invalid argument 'eq', %1", eq;
		end if;
		
		return zip(`=`, op(map(convert, leq, 'list')));
	end proc;
	protect('VFEquationToList');
	
	VF := proc (
		F :: Or(name, indexed, function),
		indices :: symbol := 'default',
		{coords :: Or(name, indexed(name)) := VectorCalculus:-GetCoordinates()},
	$)::Vector[column];
		local
			_F :: RCoefficientParameters,
			_coords :: indexed(name),
			_indices :: symbol,
			_result :: list;

		_coords := FormatCoordinates(coords);
		_F := _ParseCoefficientParameters(F);

		#` Processing indexes by default `
		_indices := indices;
		
		if ormap(verify, ['t', 'time'], _indices, 'symbol') then
			if evalb(_F:-HasArguments) then
				_F:-Arguments := [op(_F:-Arguments), 't'];
			else
				_F:-Arguments := ['__', 't'];
			end if;
			
			_indices := 'default';
		end if;
		
		if ormap(verify, ['d', 'default'], _indices, 'symbol') then
			_indices := `if`(nops(_F:-Indices) > 0, 'match', 'append');
		else
			_indices := indices;
		end if;
		
		#` Replace pattern symbol in the argument and indices lists to coordinats `
		_F:-Arguments := _ApplyTemplateCoordinates(_F:-Arguments, _coords);
		_F:-Indices := _ApplyTemplateCoordinates(_F:-Indices, _coords);

		if ormap(verify, ['n', 'none'], _indices, 'symbol') then
			_result := [_FormatCoefficientParameters(_F, _coords) $ nops(_coords)];
		elif ormap(verify, ['a', 'append'], _indices, 'symbol') then
			_result := map(
				coord -> _FormatCoefficientParameters(_F, _coords,
					'Indices' = ((l::list)->
						`if`(_F:-HasIndices and (nops(_F:-Indices) = 0),
							[], [op(l), coord])),
					'Arguments' = ((l::list)->_ApplyTemplateCoordinate(l, coord))),
				[op(_coords)]);
		elif ormap(verify, ['s', 'subs'], _indices, 'symbol') then
			if not (has('_', _F:-Indices) or has('_', _F:-Arguments)) then
				error "invalid pattern";
			end if;
			_result := map(
				coord -> _FormatCoefficientParameters(_F, _coords,
					'Indices' = ((l::list)->_ApplyTemplateCoordinate(l, coord)),
					'Arguments' = ((l::list)->_ApplyTemplateCoordinate(l, coord))),
				[op(_coords)]);
		elif ormap(verify, ['m', 'match'], _indices, 'symbol') then
			if not type(_F:-Indices, 'list(name)') or 
				nops({op(_F:-Indices)} minus {op(_coords)}) > 0 then
				error "invalid list of coordinates, %1", _F:-Indices;
			end if;

			_result := map(
				coord ->`if`(coord in {op(_F:-Indices)},
					_FormatCoefficientParameters(_F, _coords,
					'Indices' = ((l::list)->[coord]),
					'Arguments' = ((l::list)->_ApplyTemplateCoordinate(l, coord))),
				0), [op(_coords)]);
		else
			error "invalid indices: %1", indices;
		end if;
		
		return VectorCalculus:-VectorField(convert(
			_result, 'Vector'), _coords); 
	end proc;
	protect('VF');

	SubsVF := proc(
		s :: {set(equation(Vector[column])),
			list(equation(Vector[column])),
			equation(Vector[column])},
		expr :: anything,
	$)::anything;
		local
			_si :: equation(Vector[column]),
			_expr::anything;

		if type(s, Or(set, list)) then
			_expr := expr;
			for _si in s do
				_expr := SubsVF(_si, _expr);
			end do;
		else
			_expr := subs(VFEquationToList(s), expr);
		end if;
		
		return _expr;
	end proc;
	protect('SubsVF');
	
	_ParseCoefficientParameters := proc(
		F::Or(name, indexed, function),
	$)::RCoefficientParameters;
		local
			_name :: algebraic,
			_hasIndices :: boolean,
			_indices :: list,
			_hasArguments :: boolean,
			_arguments :: list;

		_hasArguments := type(F, 'function');
		if evalb(_hasArguments) then
			_name := op(0, F);
			_arguments := [op(F)];
		else
			_name := F;
			_arguments := [];
		end if;

		_hasIndices := type(_name, 'indexed');
		if evalb(_hasIndices) then
			_indices := [op(_name)];
			_name := op(0, _name);
		else
			_indices := [];
		end if;

		if not type(_name, 'name') then
			error "invalid F: %1", F;
		end if;
	
		return Record(
			':-Name' = _name,
			':-HasIndices' = _hasIndices,
			':-Indices' = _indices,
			':-HasArguments' = _hasArguments,
			':-Arguments' = _arguments);
	end proc;
	protect('_ParseCoefficientParameters');

	_FormatCoefficientParameters := proc(
		F :: RCoefficientParameters,
		coords :: indexed(name),
		{Indices :: procedure(list) := ((l::list) -> l)},
		{Arguments :: procedure(list) := ((l::list) -> l)},
	$)::Or(name, indexed, function);
		local
			_F :: RCoefficientParameters,
			_result :: Or(name, indexed, function);
		
		_F := apply(proc()::RCoefficientParameters;
			local
				_indices :: list,
				_arguments :: list;

			_indices := eval(Indices(F:-Indices));
			_arguments :=  eval(Arguments(
				`if`(nops(F:-Arguments) > 0,
					F:-Arguments,
					`if`(F:-HasArguments, [], [op(coords)]))));
			return Record(
				':-Name' = eval(F:-Name),
				':-HasIndices' = eval(nops(_indices) > 0),
				':-Indices' = eval(_indices),
				':-HasArguments' = eval(nops(_arguments) > 0),
				':-Arguments' = eval(_arguments));
		end proc);
		
		_result := _F:-Name;
		if _F:-HasIndices then
			_result := _result[op(_F:-Indices)];
		end if;

		if _F:-HasArguments then
			_result := _result(op(_F:-Arguments));
		end if;
		
		return _result;
	end proc;
	protect('_FormatCoefficientParameters');
	
	_ApplyTemplate := proc(
		items :: list,
		template :: symbol,
		replace :: anything,
	$) :: list;
		return eval(map(
			(item, t, r)->`if`(verify(item, t, 'symbol'), r, item),
			items, template, replace));
	end proc;
	protect('_ApplyTemplate');

	_ApplyTemplateCoordinate := proc(
		items :: list,
		coord :: name,
	$) :: list;
		return _ApplyTemplate(items, '_',	coord);
	end proc;
	protect('_ApplyTemplateCoordinate');
	
	_ApplyTemplateCoordinates := proc(
		items :: list,
		coords :: indexed(name),
	$) :: list;
		return _ApplyTemplate(items, '__',	'op'(coords));
	end proc;
	protect('_ApplyTemplateCoordinates');

	ModuleLoad := proc($);
		if not TypeTools:-Exists('RCoefficientParameters') then
			TypeTools:-AddType('RCoefficientParameters', 'record'(
		 		':-Name' :: name,
		 		':-HasIndices' :: boolean,
				':-Indices' :: list,
				':-HasArguments' :: boolean,
				':-Arguments'::list));
			protect('`type/RParameters`');
		end if;
	end proc;
	ModuleLoad();
	protect('ModuleLoad');

	ModuleUnload := proc($);
		TypeTools:-RemoveType('RCoefficientParameters');
	end proc;
	protect('ModuleUnload');
end module:

savelib(
	'VectorCalculusTools',
	cat(kernelopts('mapledir'), kernelopts('dirsep'), 
		"lib", kernelopts('dirsep'), 
		"VectorCalculusTools.mla")
);
