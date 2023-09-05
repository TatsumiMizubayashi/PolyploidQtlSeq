using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEffコンフィグファイル
    /// </summary>
    public class SnpEffConfigFile
    {
        /// <summary>
        /// SnpEffConfigFileの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const string DEFAULT = "";

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "sc";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "snpEffConfig";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "snpEff.config file. Not required if snpEff default config file is used.";

        /// <summary>
        /// SnpEffConfigFileを指定する。
        /// 既定のSnpEffConfigファイルを使用する場合は空文字を指定する。
        /// </summary>
        /// <param name="filePath">snpEffConfigファイルのPath</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public SnpEffConfigFile(string filePath, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Path = OptionValue.GetValue(LONG_NAME, filePath, parameterDictionary, userOptionDictionary);
            HasFile = !string.IsNullOrEmpty(Path);
        }

        /// <summary>
        /// SnpEffConfigFileを指定する。
        /// 既定のSnpEffConfigファイルを使用する場合は空文字を指定する。
        /// </summary>
        /// <param name="filePath">snpEffConfigファイルのPath</param>
        public SnpEffConfigFile(string filePath)
        {
            HasFile = !string.IsNullOrEmpty(filePath);

            if(HasFile)
            {
                if (!File.Exists(filePath)) throw new FileNotFoundException(filePath);
                Path = System.IO.Path.GetFullPath(filePath);
            }
            else
            {
                Path = "";
            }
        }

        /// <summary>
        /// SnpEffConfigファイルのPathを指定する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// SnpEffコンフィグファイルが指定されているかどうか。
        /// </summary>
        internal bool HasFile { get; }

        /// <summary>
        /// SnpEffコマンドの引数に変換する。
        /// </summary>
        /// <returns>SnpEff引数</returns>
        internal string ToSnpEffArg()
        {
            return HasFile
                ? $"-c {Path}"
                : "";
        }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Path}";
        }

    }
}
