/*
*	construct designer Javascript API. libfrag.js must also be included
* 
*	Functions: (cb = callback)
* 		cd_add_fragment(cb, cid, fid, position = last)
* 
* 		cd_rm_fragment(cb, cid, fid)
* 
* 		cd_reorder_fragment(cb, cid, [fids], [directions] = unchanged)
* 
* 		cd_get_info(cb, cid)
* 			{name, desc, length, cfs=Array[construct fragment]}
* 
* 		cd_set_info(cb, cid, name, desc)
* 
*/
/*
 * construct fragment:
 * 		- fid
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
	var d = JSON.stringify({'fid': fid});
	make_request(CD_BASE_URL + cid + '/addFragment/', d, 'Adding fragment "%s"' % fid, cb);
}

var cd_rm_fragment = function(cb, cid, fid)
{
	var d = JSON.stringify({'fid': fid});
	make_request(CD_BASE_URL + cid + '/rmFragment/', d, 'Removing fragment "%s"' % fid, cb);
}

var cd_reorder_fragments = function(cb, cid, fids, directions)
{
	if(directions != undefined)
		var d = JSON.stringify({'fids[]': fids, 'directions[]': directions});
	else
		var d = JSON.stringify({'fids[]': fids,});
		
	make_request(CD_BASE_URL + cid + '/saveOrder/', d, '', cb);
}

var cd_get_info = function(cb, cid)
{
	make_request(CD_BASE_URL + cid + '/getInfo/', undefined, 'Getting info on construct "%s"' % cid, cb)
}

var cd_set_info = function(cb, cid, name, desc)
{
	var d = JSON.stringify({'name': name, 'desc': desc,});
	make_request(CD_BASE_URL + cid + '/saveMeta/', d, 'Saving info on construct "%s"' % fid, cb);
}
