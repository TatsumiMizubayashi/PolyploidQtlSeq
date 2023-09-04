using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.Share
{
    /// <summary>
    /// リファレンスシークエンス
    /// </summary>
    public class ReferenceSequence
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "r";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "refSeq";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Reference sequence file.";

        /// <summary>
        /// リファレンスシークエンスファイルを作成する。
        /// </summary>
        /// <param name="refSeqFilePath">リファレンスシークエンスファイルPath</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public ReferenceSequence(string refSeqFilePath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, refSeqFilePath, parameterDictionary, userOptionDictionary);
            if (string.IsNullOrWhiteSpace(Path))
                throw new ArgumentException("Specify the reference sequence in the -r option.");
        }

        /// <summary>
        /// リファレンスシークエンスファイルを作成する。
        /// </summary>
        /// <param name="refSeqFilePath">リファレンスシークエンスファイルPath</param>
        public ReferenceSequence(string refSeqFilePath)
        {
            if (string.IsNullOrEmpty(refSeqFilePath)) throw new ArgumentException(null, nameof(refSeqFilePath));
            if (!File.Exists(refSeqFilePath)) throw new FileNotFoundException($"{refSeqFilePath} not found.");

            Path = refSeqFilePath;
        }



         /// <summary>
        /// リファレンスシークエンスファイルPathを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータファイルテキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }


        /// <summary>
        /// bcftools用引数に変換する。
        /// </summary>
        /// <returns>bcf引数</returns>
        internal string ToBcftoolsArg()
        {
            return $"-f {Path}";
        }
    }
}
