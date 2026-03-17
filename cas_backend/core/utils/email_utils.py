import logging
import typing as t

import sendgrid
from sendgrid.helpers.mail import Mail

from casp.services import settings

logger = logging.getLogger(__name__)

TEMPLATE_DIR = f"{settings.APP_ROOT}/casp/services/admin/templates"


class EmailSender:

    def __init__(
        self,
        from_address: str,
        sendgrid_key: t.Optional[str] = None,
        template_dir: str = TEMPLATE_DIR,
    ) -> None:
        """
        Initialize EmailSender class.

        :param sendgrid_client: SendGrid API key or client object. If this value is empty, the class will not send emails.
        :param from_address: Email address to send emails from.
        """
        if sendgrid_key is None:
            logger.warning("No SendGrid key provided.  Emails will not be sent.")
            self.sg = None
        else:
            self.sg = sendgrid.SendGridAPIClient(sendgrid_key)

        self.from_address = from_address
        self.template_dir = template_dir

    def _send_email(
        self, to: str, subject: str, content_path_html: str, content_path_plain: str, sub_values: dict[str, object] = {}
    ) -> None:
        """
        Send an email with the given content to the specified recipient.

        :param to: The email address to send the email to.
        :param subject: The subject of the email.
        :param content_path_html: The relative path to the HTML content of the email where the working directory's default is
                                ``src/casp/services/admin/templates/email``
        :param content_path_plain: The relative path to the plain text content of the email where the working directory's default is
                                ``src/casp/services/admin/templates/email``
        :param sub_values: The values to substitute into the email template.
        """
        if not self.sg:
            logger.warning("No SendGrid key provided.  Email not sent.")
            return

        with open(f"{self.template_dir}/email/{content_path_html}", "r") as f:
            html_content = f.read().format(**sub_values)
        with open(f"{self.template_dir}/email/{content_path_plain}", "r") as f:
            plain_text_content = f.read().format(**sub_values)

        message = Mail(
            from_email=self.from_address,
            to_emails=to,
            subject=subject,
            html_content=html_content,
            plain_text_content=plain_text_content,
        )
        self.sg.send(message)

    def send_welcome_email(self, email: str, key: str) -> None:
        """
        Send a welcome email to the given email address with the given key.

        :param email: The email address to send the email to.
        :param key: The key to include in the email.
        """
        self._send_email(
            to=email,
            subject="Welcome to the Cellarium Cell Annotation Service",
            content_path_html="welcome.html",
            content_path_plain="welcome.txt",
            sub_values={"key": key},
        )

    def send_new_key_email(self, email: str, key: str) -> None:
        """
        Send an email to the given email address with the given key.

        :param email: The email address to send the email to.
        :param key: The key to include in the email.
        """
        self._send_email(
            to=email,
            subject="New Cellarium Cell Annotation Service Key you Requested",
            content_path_html="new_key.html",
            content_path_plain="new_key.txt",
            sub_values={"key": key},
        )
