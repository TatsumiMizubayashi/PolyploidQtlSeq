namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Mappingサンプル設定
    /// </summary>
    internal class MappingSampleSettings
    {
        /// <summary>
        /// Mappngオプションを作成する。
        /// </summary>
        /// <param name="settingValue">Mappingサンプル設定値</param>
        public MappingSampleSettings(IMappingSampleSettingValue settingValue)
        {
            Parent1Directory = new Parent1Directory(settingValue.Parent1Dir);
            Parent2Directory = new Parent2Directory(settingValue.Parent2Dir);
            Bulk1Directory = new Bulk1Directory(settingValue.Bulk1Dir);
            Bulk2Directory = new Bulk2Directory(settingValue.Bulk2Dir);
        }


        /// <summary>
        /// Parent1ディレクトリ
        /// </summary>
        public Parent1Directory Parent1Directory { get; }

        /// <summary>
        /// Parent2ディレクトリ
        /// </summary>
        public Parent2Directory Parent2Directory { get; }

        /// <summary>
        /// Bulk1ディレクトリ
        /// </summary>
        public Bulk1Directory Bulk1Directory { get; }

        /// <summary>
        /// Bulk2ディレクトリ
        /// </summary>
        public Bulk2Directory Bulk2Directory { get; }
    }
}
