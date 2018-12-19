function [decoding_ind,decoding_start,decoding_stop] = get_decoding_ind(ts,num_decode_chunks,num_int)

dt = median(diff(ts))/1e6;
buffer = ceil(1*60/dt); % time bins towards end or start of session that I toss out
decoding_start = round(linspace(buffer,numel(ts) - buffer,num_decode_chunks));
decoding_stop = decoding_start+num_int;

decoding_ind = [];
for k = 1:num_decode_chunks
    decode_chunk_k = (decoding_start(k):decoding_start(k)+num_int)';
    decoding_ind = [decoding_ind; decode_chunk_k];
end

return
