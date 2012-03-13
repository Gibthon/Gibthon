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
* 		cd_get_info(cb, cid, ecb)
* 			cb = function({name, desc, Array[fragment] fs, Array[construct fragment] cfs})
* 
* 		cd_set_info(cb, cid, name, desc, ecb)
* 			cb = function()
* 
*/
/*
 * construct fragment:
 * 		d = json data
 * 		frag = fragment, fetched from server if undefined or wrong
 * 
 * */

function ConstructFragment(d, f)
{	
	if(d == undefined) d = {};
	
	this.id = d.cfid;
	this.strand = d.direction;
	this.s_offset = d.s_offset;
	this.s_feat = d.s_feat;
	this.e_offset = d.e_offset;
	this.e_feat = d.e_feat;
	this.order = d.order;
	
	if(((f == undefined) || (f.fid != d.fid)) && (d.fid != undefined))
	{
		get_meta(d.fid, function(d) {
			f = new Fragment(d);
		});
	}
	this.isValid = function() {return f != undefined;};
	this.startPos = function()
	{
		if(this.id == undefined) return 0;
		//console.log('startPos: s_feat: '+this.s_feat+' s_offset: '+this.s_offset);
		if(this.s_feat > 0)
		{
			var sf = f.getFeatById(this.s_feat);
			if(sf != null)
				return this.s_offset + sf.start;
		}
		return this.s_offset;
	};
	this.endPos = function()
	{
		if(this.id == undefined) return f.length;
		if(this.e_feat > 0)
		{
			var ef = f.getFeatById(this.e_feat);
			if(ef != null)
				return this.e_offset + ef.end;
		}
		return f.length - this.s_offset;
	};
	this.length = function()
	{
		return Math.abs(this.endPos() - this.startPos());
	};
}

var CD_BASE_URL = '/gibthon/api/'

var cd_add_fragment = function(cb, cid, fid, position, ecb)
{
	make_request(	CD_BASE_URL + cid + '/addFragment/', 
					JSON.stringify({'fid': fid}), 
					'Adding fragment "' + fid + '"', 
					cb);
}

var cd_rm_fragment = function(cb, cid, cfid, ecb)
{
	make_request(	CD_BASE_URL + cid + '/rmFragment/', 
					JSON.stringify({'cfid': cfid,}), 
					'Removing ConstructFragment "' + cfid + '"', 
					cb);
}

var cd_reorder_fragments = function(cb, cid, d, ecb)
{
	make_request(	CD_BASE_URL + cid + '/saveOrder/', 
					JSON.stringify({'d[]': d,}), 
					'Reordering fragments', 
					cb);
}

var cd_get_info = function(cb, cid, ecb)
{
	make_request(	CD_BASE_URL + cid + '/getInfo/', 
					undefined, 
					'Getting info on construct "'+cid+'"', 
					cb);
}

var cd_set_info = function(cb, cid, name, desc, ecb)
{
	make_request(	CD_BASE_URL + cid + '/saveMeta/', 
					JSON.stringify({'name': name, 'desc': desc,}), 
					'Saving info on construct "'+fid+'"', 
					cb);
}
