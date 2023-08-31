using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 入力RawFastqディレクトリ
    /// </summary>
    public class InputRawFastqDirectory
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
        public const string LONG_NAME = "inputDir";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Raw fastq directory.";

        /// <summary>
        /// 入力RawFastqディレクトリを作成する。
        /// </summary>
        /// <param name="inputDirPath">入力RawFastqディレクトリPath</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        [Obsolete("削除予定")]
        public InputRawFastqDirectory(string inputDirPath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, inputDirPath, parameterDictionary, userOptionDictionary);

            if (string.IsNullOrWhiteSpace(Path)) throw new ArgumentException("Specify the Raw Fastq directory in the -i option.");
            if (!Directory.Exists(Path)) throw new DirectoryNotFoundException(Path);
        }

        /// <summary>
        /// 入力RawFastqディレクトリのPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        [Obsolete("削除予定")]
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }

        /// <summary>
        /// ディレクトリ内のFastqファイルをFastqFilePair情報に変換する。
        /// </summary>
        /// <returns>Fastqファイルペア配列</returns>
        internal FastqFilePair[] ToFastqFilePairs()
        {
            return FastqFilePairEnumerator.Enumerate(Path);
        }
    }
}
