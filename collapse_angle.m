% COLLAPSE_ANGLE_IQ
%
% IQ = collapse_angle_IQ( IQ_raw, num_angles )
%
% This function processes angle-compounded IQ data.

function IQ = collapse_angle( IQ_raw, num_angles )
	IQ = nan( size(IQ_raw) - [0 0 num_angles-1] );
	for frame = 1:size(IQ,3)
		IQ(:,:,frame) = mean( IQ_raw(:,:,frame:frame+num_angles-1), 3 );
	end
end
