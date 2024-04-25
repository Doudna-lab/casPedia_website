# == Native Modules

# == Installed Modules
import os
from google.oauth2 import service_account
from googleapiclient.discovery import build
import chardet
# == Project Modules


def copy_doc_with_comments(service, document_id):
  """
  Creates a copy of a Google Doc with comments included.

  Args:
	  service: The authorized Docs API service object.
	  document_id: The ID of the template document to copy.

  Returns:
	  The ID of the newly created copy.
  """
  try:
	# Set copy behavior to include comments and suggestions
	body = {
	  'copyCommentsAndSuggestions': True
	}

	# Create a new copy of the document
	copied_doc = service.documents().copy(body=body, documentId=document_id).execute()
	return copied_doc.get('id')

  except HttpError as error:
	print(f"An error occurred: {error}")
	return None


def create_copy_of_template(drive_service, sheets_service, parent_folder_id, template_file_id, copy_name):
	# Retrieve content of template Google Docs file
	# Retrieve content of template Google Sheets file
	template_response = drive_service.files().get(fileId=template_file_id, fields="*").execute()

	# Create a new Google Docs file
	new_file_metadata = {
		'name': copy_name,
		'mimeType': 'application/vnd.google-apps.spreadsheet'
	}

	# Set the parent folder if provided
	if parent_folder_id:
		new_file_metadata['parents'] = [parent_folder_id]

	new_file = drive_service.files().create(body=new_file_metadata).execute()
	new_sheet_id = new_file['id']

	body = {
		'requests': [{
			'copyPaste': {
				'source': {
					'sheetId': template_file_id  # Provide the source sheet ID
				},
				'destination': {
					'sheetId': 0  # Specify the destination sheet ID, e.g., 0 for the first sheet
				},
				'pasteType': 'PASTE_NORMAL'  # Specify the paste type, e.g., 'PASTE_NORMAL'
			}
		}]
	}

	sheets_service.spreadsheets().batchUpdate(spreadsheetId=new_sheet_id, body=body).execute()

	return template_response


# AIzaSyAiJetZN02EtxGxTCjv1cGi5lkfpwO6-Zg
def main():
	# Authentication
	SCOPES = ['https://www.googleapis.com/auth/drive',
			  'https://www.googleapis.com/auth/documents',
			  'https://www.googleapis.com/auth/spreadsheets']
	SERVICE_ACCOUNT_FILE = '/home/ubuntu/casPedia_website/keys/service_account.json'

	credentials = service_account.Credentials.from_service_account_file(
		SERVICE_ACCOUNT_FILE, scopes=SCOPES)

	# Create a Google Drive and Google Docs service
	drive_service = build('drive', 'v3', credentials=credentials)
	docs_service = build('docs', 'v1', credentials=credentials)
	sheets_service = build('sheets', 'v4', credentials=credentials)

	template_file_id = '1AO0v1u2vh3sWTgwDDifd0ixRW6PyZwfZgrIhtOa9rqU'  # Provide the ID of your template Google Docs file
	parent_folder_id = '1WrOK8D270zQ4lSGpg0X6rHtacDOrANG0'
	copy_name = 'COPIED_WITH_API'  # Provide the desired name for the copied file
	new_file_id = create_copy_of_template(drive_service, sheets_service, parent_folder_id, template_file_id, copy_name)
	print(f"New file created with ID: {new_file_id}")


if __name__ == "__main__":
	main()
