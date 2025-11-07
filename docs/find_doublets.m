% Look for shots in shotlist and return a structure with:
% Run: date
% Shot: shot number
% Topic: we are only selecting PDJ entries, they usually log blips
% Text: the actual log (capped at 5000 characters)

START_SHOT = 87700; % May 2025
END_SHOT = 87899; % Sept 2025

shotlist = START_SHOT:END_SHOT; % for small batches of shot u can do a single call
% topic_list = "('PdJ')";
topic_list = "('SL')";

doublets = [];

for shot = shotlist
	ddd = sprintf('SELECT Run, Shot, Username, Topic, Text, Entered FROM Entries WHERE (Shot in (%d%s)) AND (Topic in %s) and (voided is NULL) ORDER BY Shot DESC, Topic limit 5000',shot,sprintf(',%d',shot),topic_list);
	F = sqllogbook(ddd);
	try
		% print detected doublet numbers
		% fprintf('Text: %s \n', F.Text{:})
		for i = 1 : numel(F.Text)
			if contains(F.Text(i),'doublet', 'IgnoreCase', true)
				fprintf('Shot %d -> doublet\n',F.Shot(i))
				doublets = [doublets, shot];
			end
		end
	catch ME
		fprintf('Error processing shot %d: %s \n', shot, ME.message);
        % for k = 1:length(ME.stack)
        %     fprintf('\t\tIn %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
        % end
	end
end

fprintf('Doublet shots detected: %s\n', mat2str(doublets));

















