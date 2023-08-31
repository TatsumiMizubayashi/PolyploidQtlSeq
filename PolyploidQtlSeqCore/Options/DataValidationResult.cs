namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// データ検証結果
    /// </summary>
    public class DataValidationResult
    {
        private static readonly string _separator = new('-', 60);

        /// <summary>
        /// エラーなしデータ検証結果インスタンスを作成する。
        /// </summary>
        public DataValidationResult()
        {
            HasError = false;
            OptionName = "";
            ErrorMessage = "";
        }

        /// <summary>
        /// エラー情報を持つデータ検証結果インスタンスを作成する。
        /// </summary>
        /// <param name="shortName">Option ShortName</param>
        /// <param name="longName">Option LongName</param>
        /// <param name="errorMessage">エラーメッセージ</param>
        public DataValidationResult(string shortName, string longName, string errorMessage)
        {
            HasError = true;
            OptionName = $"-{shortName}|--{longName}";
            ErrorMessage = errorMessage;
        }

        /// <summary>
        /// エラーの有無を取得する。
        /// </summary>
        public bool HasError { get; }

        /// <summary>
        /// オプションスイッチのShort/LongNameを取得する。
        /// </summary>
        public string OptionName { get; }

        /// <summary>
        /// エラーメッセージを取得する。
        /// </summary>
        public string ErrorMessage { get; }

        /// <summary>
        /// エラーメッセージを標準出力に出力する。
        /// </summary>
        public void Print()
        {
            if (!HasError) return;

            Console.Error.WriteLine($"Option: {OptionName}");
            Console.Error.WriteLine(ErrorMessage);
            Console.Error.WriteLine(_separator);
        }
    }
}
