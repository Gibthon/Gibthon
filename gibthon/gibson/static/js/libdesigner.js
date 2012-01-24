/*
*	construct designer Javascript API. libfrag.js must also be included
* 
*	Functions: (cb = callback)
* 		cd_add_fragment(cb, cid, fid, position = last) 
* 			cb = function(construct fragment)
* 
* 		cd_rm_fragment(cb, cid, cfid)
* 			cb = function()
* 
* 		cd_reorder_fragment(cb, cid, {Array[] cfid, Array[] direction})
* 			cb = function()
* 
* 		cd_get_info(cb, cid)
* 			cb = function({name, desc, length, Array[construct fragment] cfs})
* 
* 		cd_set_info(cb, cid, name, desc)
* 			cb = function()
* 
*/
/*
 * construct fragment:
 * 		- cfid
 * 		- fid
 * 		- order
 * 		- direction [1,-1]
 * 		- s_feat: start feature, null or feature ID
 * 		- s_offset: start offset, integer
 * 		- e_feat: end feature
 * 		- e_offset: end
 * 
 * */
 
var CD_BASE_URL = '/gibthon/api/'

var cd_add_fragment = function(cb, cid, fid, position)
{
	make_request(	CD_BASE_URL + cid + '/addFragment/', 
					JSON.stringify({'fid': fid}), 
					'Adding fragment "' + fid + '"', 
					cb);
}

var cd_rm_fragment = function(cb, cid, cfid)
{
	make_request(	CD_BASE_URL + cid + '/rmFragment/', 
					JSON.stringify({'cfid': cfid}), 
					'Removing ConstrictFragment "' + cfid + '"', 
					cb);
}

var cd_reorder_fragments = function(cb, cid, d)
{
	make_request(	CD_BASE_URL + cid + '/saveOrder/', 
					JSON.stringify({'d[]': d,}), 
					'Reordering fragments', 
					cb);
}

var cd_get_info = function(cb, cid)
{
	make_request(	CD_BASE_URL + cid + '/getInfo/', 
					undefined, 
					'Getting info on construct "'+cid+'"', 
					cb);
}

var cd_set_info = function(cb, cid, name, desc)
{
	make_request(	CD_BASE_URL + cid + '/saveMeta/', 
					JSON.stringify({'name': name, 'desc': desc,}), 
					'Saving info on construct "'+fid+'"', 
					cb);
}
