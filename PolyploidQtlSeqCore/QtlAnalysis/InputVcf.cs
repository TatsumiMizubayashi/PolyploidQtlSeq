using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 入力VCFファイル
    /// </summary>
    public class InputVcf
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "i";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "inputVcf";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Input VCF file.";

        /// <summary>
        /// 入力ファイルを作成する。
        /// </summary>
        /// <param name="filePath">VCFファイルPath</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public InputVcf(string filePath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, filePath, parameterDictionary, userOptionDictionary);

            if (string.IsNullOrWhiteSpace(Path)) throw new ArgumentException("Specify a VCF file for the -i option.");
        }

        /// <summary>
        /// 入力ファイルを作成する。
        /// </summary>
        /// <param name="filePath">VCFファイルPath</param>
        public InputVcf(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentNullException(nameof(filePath));

            Path = System.IO.Path.GetFullPath(filePath);
        }

        /// <summary>
        /// VCFファイルパスを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns></returns>
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }

    }
}
